import sys
import os
import yaml
import astropy.table as at
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from SNmodel import spline_utils


mask_bit_locations = {'verbose':1,'dump':2}
DEFAULT_BAYESN_MODEL='M20'
ALLOWED_BAYESN_MODEL=['M20', 'T21']
PRODUCTS_DIR = os.getenv('PRODUCTS')
BAYESN_MODEL_DIR = os.path.join(PRODUCTS_DIR, 'bayesn', 'SNmodel', 'model_files')
BAYESN_MODEL_COMPONENTS = ['l_knots', 'L_Sigma_epsilon', 'M0_sigma0_RV_tauA', 'tau_knots', 'W0', 'W1']


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

            print('PARAMS FILE:')
            print(self.params_file_contents)

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
            self.parameter_names = ['THETA1','AV','RV','DELTAM','EPSILON','TMAX']
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


    def fetchSED_BAYESN(self,trest,maxlam=5000,external_id=1,new_event=1,hostpars=''):
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
            self.parameter_values = {key:0.+0.1*external_id for key in self.parameter_names}

		#Matrix needed for wavelength spline
        # KD_l = spline_utils.invKD_irr(self.l_k)

        ind =  np.abs(self.phase - trest).argmin()
        ind_flux = ind*self.wavelen
        flux = self.flux[ind_flux:ind_flux+self.wavelen]
        return list(flux)



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
  mySED=gensed_BAYESN('$WFIRST_USERS/jpierel/pySED_test/SNEMO.P20/',2,[],'z,AGE,ZCMB,METALLICITY')


if __name__=='__main__':
    main()
