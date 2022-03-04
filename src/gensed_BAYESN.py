import sys
import os
import yaml
import astropy.table as at
import extinction
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from SNmodel import spline_utils


mask_bit_locations = {'verbose':1,'dump':2}
DEFAULT_BAYESN_MODEL='M20'
ALLOWED_BAYESN_MODEL=['M20', 'T21']
PRODUCTS_DIR = os.getenv('PRODUCTS')
BAYESN_MODEL_DIR = os.path.join(PRODUCTS_DIR, 'bayesn', 'SNmodel', 'model_files')
BAYESN_MODEL_COMPONENTS = ['l_knots', 'L_Sigma_epsilon', 'M0_sigma0_RV_tauA', 'tau_knots', 'W0', 'W1']
#ST: Computes a stupid nuisance factor
GAMMA = np.log(10)/2.5


def print_err():
  print(".·´¯`(>▂<)´¯`·.  ┻━┻︵ \(°□°)/ ︵ ┻━┻ ABORT Python on Fatal Error.")
  raise RuntimeError


class gensed_BAYESN:
    def __init__(self,PATH_VERSION,OPTMASK,ARGLIST,HOST_PARAM_NAMES):

        try:
            self.verbose = OPTMASK & (1 << mask_bit_locations['verbose']) > 0
            self.host_param_names = [x.upper() for x in HOST_PARAM_NAMES.split(',')]
            self.dump = OPTMASK & (1 << mask_bit_locations['dump'])>0
            self.PATH_VERSION = os.path.expandvars(PATH_VERSION)

            # check if a param file exists
            self.paramfile = None
            param_files = ['bayesn.params','BAYESN.PARAMS', 'BAYESN.params']
            for param_file in param_files:
                if os.path.exists(os.path.join(self.PATH_VERSION, param_file)):
                    self.paramfile = os.path.join(self.PATH_VERSION,param_file)
                    break
            if self.paramfile is None:
                raise RuntimeError(f'param file not found! in {self.PATH_VERSION}. Looking for one of {param_files}')

            self.params_file_contents = yaml.load(open(self.paramfile),
                                                  Loader=yaml.FullLoader)

            SNDATA_PATH = os.getenv('SNDATA_ROOT')
            if SNDATA_PATH is None:
                raise RuntimeError('SNDATA_ROOT is not defined! Check env!')
            hsiao_model = os.path.join(SNDATA_PATH,'snsed','Hsiao07.dat')
            if not os.path.isfile(hsiao_model):
                raise RuntimeError(f'Cannot load Hsiao Model - check if {hsiao_model} exists.')

            self._hsiao = at.Table.read(hsiao_model, format='ascii', names=('phase','wave','flux'))

            ### FILL IN THESE REQUIRED ELEMENTS
            self.phase = np.unique(self._hsiao['phase'])
            self.wave = np.unique(self._hsiao['wave'])
            self.wavelen = len(self.wave)
            self.flux = self._hsiao['flux']
            self.parameter_names = ['THETA1','AV','RV','DELTAM','EPSILON','TMAX', 'REDSHIFT']
            self.parameter_values = {key:-9. for key in self.parameter_names}

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print('Python Error :',e)
            print('gensed_BAYESN.py, line number: %i'%exc_tb.tb_lineno)
            print_err()

        # load the bayesn model files
        BAYESN_MODEL = self.params_file_contents.get('BAYESN_MODEL', DEFAULT_BAYESN_MODEL)
        if BAYESN_MODEL not in ALLOWED_BAYESN_MODEL:
            print(f'gensed_BAYESN.py: INVALID BAYESN_MODEL {BAYESN_MODEL} - must be one of {ALLOWED_BAYESN_MODEL}')
            print_err()

        bayesn_model_dir = os.path.join(BAYESN_MODEL_DIR, f'{BAYESN_MODEL}_model')
        self._bayesn_components = {comp:np.genfromtxt(os.path.join(bayesn_model_dir, f'{comp}.txt')) for comp in BAYESN_MODEL_COMPONENTS}

        #ST: Computes spline invrse KD matrices.
        self.KD_t = spline_utils.invKD_irr(self._bayesn_components["tau_knots"])
        self.KD_l = spline_utils.invKD_irr(self._bayesn_components["l_knots"])
        self.J_l = spline_utils.spline_coeffs_irr(self.wave, self._bayesn_components["l_knots"], self.KD_l)

        #ST: Extracts the M0 parameter (this is kind of horrible)
        self.M0 = self._bayesn_components["M0_sigma0_RV_tauA"][0]

        return


    def fetchSED_NLAM(self):
        """
        Returns the length of the wavelength vector
        """
        return self.wavelen


    def fetchSED_LAM(self):
        """
        Returns the wavelength vector
        """
        return list(self.wave)


    def fetchSED_BAYESN(self, trest, maxlam=5000, external_id=1, new_event=1, hostpars=''):
        """
        Returns the flux at every wavelength, for a given phase.

        Parameters
        ----------
        trest : float
             The rest frame phase at which to calculate the flux
        maxlam : int
             A maximum number of wavelength bins. If your wavelength
             vector is longer than this number (default is arbitrary),
             program should abort
        external_id : int
             ID for SN
        new_event : int
             1 if new event, 0 if same SN
        hostpars : str
             Comma separated list of host parameters

        Returns
        -------
        A list of length self.wavelen containing the flux at
        every wavelength in self.wave, at the phase trest
        """
        if new_event != 1:
            newSN=False
        else:
            newSN=True
            # double or assymetric gaussian
            # genPDF map - generic multi-D map - captures correlations between parameters
            # can be any function
            # SNANA - sec 4.3.2
            # snlc_sim.c has header
            # get_random_genPDF
            print(self.host_param_names)
            print(hostpars)
            useful_pars = {x:hostpars[i] for i, x in enumerate(self.host_param_names)}
            print(useful_pars, 'kwargs')
            self.setParVals_BAYESN(**useful_pars)
            print(self.parameter_values)

            #self.parameter_values = {key:0.+0.1*external_id for key in self.parameter_names}

        #ST: Three lines originally from GSN
        ind =  np.abs(self.phase - trest).argmin()
        ind_flux = ind*self.wavelen
        flux = self.flux[ind_flux:ind_flux+self.wavelen]

        ########## NEW CODE FROM ST BELOW HERE (FOR GUIDANCE) ##########


        #ST: Computes matrices that do interpolation
        #    Probably can't be precomputed
        #    Assumes that `trest` is a float, and `self.wave` is a 1D
        #    list or numpy array of rest frame wavelengths
        J_t = spline_utils.spline_coeffs_irr([trest], self._bayesn_components["tau_knots"], self.KD_t).T

        #ST: Computes host extinction
        #    This assumes we can use the Kyle Barbary extinction.py package
        #    If we can't, I have the necessary code for this
        #    Note: This may need `self.wave` to be converted to a numpy
        #    array, if it is a list.
        R_host = extinction.fitzpatrick99(self.wave, self.parameter_values["AV"], self.parameter_values["RV"])


        #ST: Computes the spline knots
        W = self._bayesn_components["W0"] + self.parameter_values["THETA1"]*self._bayesn_components["W1"] + self.parameter_values["EPSILON"]

        #ST: Interpolates to `trest` and `self.wave`
        #    If we have done this right, this should be the same length
        #    as `flux`
        JWJ = np.linalg.multi_dot([self.J_l, W, J_t]).squeeze()


        #ST: Multiplies correction into Hsiao fluxes
        #    Stilde is essentially the host-dust-extinguished
        #    rest-frame SED, without the M0 and DELTAM normalisation
        #    factors
        S_tilde = flux*np.exp(-GAMMA*(JWJ + R_host))

        #ST: Applies the constant (in time and wavelength) factors
        flux = np.exp(-GAMMA*(self.M0 + self.parameter_values["DELTAM"]))*S_tilde

        ########## ORIGINAL RETURN STATEMENT FROM GSN ##########
        return list(flux)


    def setParVals_BAYESN(self, **kwargs):
        """
        Sets the values of parameter names
        """
        print(kwargs, 'set')
        keynames = kwargs.keys()
        for key in keynames:
            if key not in self.parameter_names:
                message = f'Unknown parameter {key}'
                print(message)
                print_err()
        self.parameter_values.update(kwargs)


    def fetchParNames_BAYESN(self):
        """
        Returns the names of model parameters
        """
        return list(self.parameter_names)


    def fetchNParNames_BAYESN(self):
        """
        Returns the number of model parameters
        """
        return len(self.parameter_names)


    def fetchParVals_BAYESN_4SNANA(self,varname):
        """
        Returns the value of parameter 'varname'

        Parameters
        ----------
        varname : str
             A parameter name from self.parameter_names
        """
        return self.parameter_values[varname]


def main():
    mySED=gensed_BAYESN('$SNANA_LSST_USERS/gnarayan/bayesn/',2,[],'z,AGE,ZCMB,METALLICITY')

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(111)
    mySED.setParVals_BAYESN(THETA1=0., AV=0.1 , RV=3.1 , DELTAM=0. , EPSILON=0. , TMAX=0.)
    flux = mySED.fetchSED_BAYESN(0, new_event=2)
    ax1.plot(mySED.wave, np.array(flux), 'g-')

    mySED.setParVals_BAYESN(THETA1=2., AV=0.3 , RV=3.1 , DELTAM=0. , EPSILON=0. , TMAX=0.)
    flux = mySED.fetchSED_BAYESN(0, new_event=2)
    ax1.plot(mySED.wave, np.array(flux), 'r-')

    mySED.setParVals_BAYESN(THETA1=-1.7, AV=0. , RV=3.1 , DELTAM=0.5 , EPSILON=0. , TMAX=0.)
    flux = mySED.fetchSED_BAYESN(0, new_event=2)
    ax1.plot(mySED.wave, np.array(flux), 'b-')

    ax1.plot(mySED.wave, np.array(flux))
    fig.savefig('test.png')

    fig2 = plt.figure(figsize=(10, 20))
    ax2 = fig2.add_subplot(111)
    ind = (mySED.wave >= 3000) & (mySED.wave <= 7000)
    for trest in np.linspace(-10, 40, 20):

        mySED.setParVals_BAYESN(THETA1=-1.7, AV=0. , RV=3.1 , DELTAM=0. , EPSILON=0. , TMAX=0.)
        flux = mySED.fetchSED_BAYESN(trest, new_event=2)
        ax2.plot(mySED.wave[ind], np.array(flux)[ind]/np.array(flux)[ind].max() + trest, 'b-')

        mySED.setParVals_BAYESN(THETA1=0., AV=0.1 , RV=3.1 , DELTAM=0. , EPSILON=0. , TMAX=0.)
        flux = mySED.fetchSED_BAYESN(trest, new_event=2)
        ax2.plot(mySED.wave[ind], np.array(flux)[ind]/np.array(flux)[ind].max() + trest, 'g-')

        mySED.setParVals_BAYESN(THETA1=2., AV=0.3 , RV=3.1 , DELTAM=0. , EPSILON=0. , TMAX=0.)
        flux = mySED.fetchSED_BAYESN(trest, new_event=2)
        ax2.plot(mySED.wave[ind], np.array(flux)[ind]/np.array(flux)[ind].max() + trest, 'r-')
    fig2.savefig('phase_sequence.png')






if __name__=='__main__':
    main()
