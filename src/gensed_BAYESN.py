import sys
import os
import yaml
import astropy.table as at


mask_bit_locations = {'verbose':1,'dump':2}


def print_err():
  print(".·´¯`(>▂<)´¯`·.  ┻━┻︵ \(°□°)/ ︵ ┻━┻ ABORT Python on Fatal Error.")
  raise RuntimeError


class gensed_BAYESN:
    def __init__(self,PATH_VERSION,OPTMASK,ARGLIST,HOST_PARAM_NAMES):
        try:
            self.verbose = OPTMASK & (1 << mask_bit_locations['verbose']) > 0
            self.host_param_names = [x.upper() for x in HOST_PARAM_NAMES.split(',')]
            self.dump = OPTMASK & (1 << mask_bit_locations['dump'])>0
            self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))

            # check if a param file exists
            if os.path.exists(os.path.join(self.PATH_VERSION,'bayesn.params')):
                self.paramfile = os.path.join(self.PATH_VERSION,'bayesn.params')
            # look for all upper case extension as well as lower case
            elif os.path.exists(os.path.join(self.PATH_VERSION,'BAYESN.PARAMS')):
                self.paramfile = os.path.join(self.PATH_VERSION,'BAYESN.PARAMS')
            else:
                raise RuntimeError('param file %s not found!'%\
                                   os.path.join(self.PATH_VERSION,'BAYESN.params'))

            #self.params_file_contents = yaml.load(open(self.paramfile),
            #                                      Loader=yaml.FullLoader)

            #print('PARAMS FILE:')
            #print(self.params_file_contents)

            SNDATA_PATH = os.getenv('SNDATA_ROOT')
            if SNDATA_PATH is None:
                raise RuntimeError('SNDATA_ROOT is not defined! Check env!')
            hsiao_model = os.path.join(SNDATA_PATH,'snsed','Hsiao07.dat')
            if not os.path.isfile(hsiao_model):
                raise RuntimeError(f'Cannot load Hsiao Model - check if {hsiao_model} exists.')

            self._hsiao = at.Table.read(hsiao_model, format='ascii', names=('phase','wave','flux'))
            print(self._hsiao)

            ### FILL IN THESE REQUIRED ELEMENTS
            #self.wave =
            #self.wavelen = len(self.wave)
            #self.parameter_names = []
            #self.parameter_values = {}

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print('Python Error :',e)
            print('gensed_BAYESN.py, line number: %i'%exc_tb.tb_lineno)
            print_err()
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


    def fetchSED_SNEMO(self,trest,maxlam=5000,external_id=1,new_event=1,hostpars=''):
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

        try:
            if len(self.wave)>maxlam:
                raise RuntimeError("Your wavelength array\
                        cannot be larger than %i but is %i"%(maxlam,len(self.wave)))

            # Calculation for flux at self.wave wavelengths for phase trest
            return flux_at_all_wavelengths_for_trest
        except Exception as e:
            print('Python Error :',e)
            print_err()


    def fetchParNames_SNEMO(self):
        """
        Returns the names of model parameters
        """
        return list(self.parameter_names)


    def fetchNParNames_SNEMO(self):
        """
        Returns the number of model parameters
        """
        return len(self.self.parameter_names)


    def fetchParVals_SNEMO_4SNANA(self,varname):
        """
        Returns the value of parameter 'varname'

        Parameters
        ----------
        varname : str
             A parameter name from self.parameter_names
        """
        return self.parameter_values[varname]


def main():
  mySED=gensed_SNEMO('$WFIRST_USERS/jpierel/pySED_test/SNEMO.P20/',2,[],'z,AGE,ZCMB,METALLICITY')


if __name__=='__main__':
    main()
