import sys
import os
# disable threading
os.environ["MKL_NUM_THREADS"]        = "1"
os.environ["NUMEXPR_NUM_THREADS"]    = "1"
os.environ["OMP_NUM_THREADS"]        = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS']   = "1"
import re
import yaml
import astropy.table as at
import extinction
import numpy as np
import scipy.stats as st
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt


mask_bit_locations = {'verbose':1,'dump':2}
DEFAULT_BAYESN_MODEL='M20'
ALLOWED_BAYESN_MODEL=['M20', 'T21']
ALLOWED_BAYESN_PARAMS=['THETA1','DELTAM'] # I suppose we could allow EPSILON as well...
ALLOWED_SNANA_DISTRIBUTION_KEYS=['PEAK','SIGMA','RANGE']
PRODUCTS_DIR = os.getenv('SNDATA_ROOT')
BAYESN_MODEL_DIR = os.path.join(PRODUCTS_DIR, 'models', 'bayesn')
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
            self.parameter_names = ['THETA1','AV','RV','DELTAM','TMAX', 'REDSHIFT']
            self.parameter_values = {key:0. for key in self.parameter_names}
            self.parameter_values['RV'] = 3.1 # make sure Rv has a sane default

            ### SETUP THE BAYESN PARAMETER DISTRIBUTIONS
            ### THESE WILL BE SAMPLED TO UPDATE parameter_values AS NEEDED
            self.parse_bayesn_param_keys()
            self.setup_bayesn_param_distributions()

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

        bayesn_model_dir = os.path.join(BAYESN_MODEL_DIR, f'BAYESN.{BAYESN_MODEL}')
        self._bayesn_components = {comp:np.genfromtxt(os.path.join(bayesn_model_dir,\
                                    f'{comp}.txt')) for comp in BAYESN_MODEL_COMPONENTS}

        self._nepsilon = len(self._bayesn_components["L_Sigma_epsilon"])

        # GN - initalize the epsilon vector
        self.parameter_names += [f'EPSILON{i:02d}' for i in range(self._nepsilon)]
        for i in range(self._nepsilon):
            self.parameter_values[f'EPSILON{i:02d}'] = 0.

        #ST: Computes spline invrse KD matrices.
        self.KD_t = invKD_irr(self._bayesn_components["tau_knots"])
        self.KD_l = invKD_irr(self._bayesn_components["l_knots"])
        self.J_l =  spline_coeffs_irr(self.wave, self._bayesn_components["l_knots"], self.KD_l)

        #ST: Extracts the M0 parameter (this is kind of horrible)
        self.M0 = self._bayesn_components["M0_sigma0_RV_tauA"][0]

        return


    def parse_bayesn_param_keys(self, arglist=None):
        """
        We allow users to override the default BayeSN parameter
        ranges by specifying the ARGLIST variable in the SNANA
        input file, but we have to parse that single string

        The keys in the string must  mirror those in bayesn.params
        (where they can mercifully be specified as separate lines)
        and must obey SNANA conventions
        i.e. the keys must be in the form:
        GEN[RANGE|PEAK|SIGMA]_<PARAM_NAME>: scalar or space-separated list

        Sets a dictionary bayesn_distrib_params with the mapping of
        key:value with value parsed as a scalar or list of floats

        and a dictionary bayesn_param_config with the mapping of
        PARAM_NAME: list of SNANA keys set (RANGE|PEAK|SIGMA)
        """

        # the pattern is YAML format - key: val pairs
        # however, unlike YAML, we get all the args in
        # a single arglist
        if arglist is not None:
            keyvals = re.split('(\w+):', arglist)[1:]
            # key vals come in pairs - even indices are keys
            # odd indices are vals
            keys = keyvals[::2]
            vals = keyvals[1::2]
            self.bayesn_distrib_params.update(dict(zip(keys, vals)))
        else:
            self.bayesn_distrib_params = {key:val for key,val in self.params_file_contents.items()
                                                if key.startswith('GEN')}

        self.bayesn_param_config = {}
        # we need to do some input validating
        for key, val in self.bayesn_distrib_params.items():

            # check the keywords are valid
            if not key.startswith('GEN'):
                message = f'Do not recognize key {key} in ARGLIST. Keys must begin with GEN'
                raise RuntimeError(message)
            keysuffix, paramname = key.replace('GEN','').split('_')

            # looks like distributions in SNANA have to be PEAK/SIGMA/RANGE/SKEW
            # I don't understand how SKEW is implemented so skipping that for now
            if keysuffix not in ALLOWED_SNANA_DISTRIBUTION_KEYS:
                message = f'Do not recognize key {key} in ARGLIST.\n'
                message += f'Keys must contain f{ALLOWED_SNANA_DISTRIBUTION_KEYS}.'
                raise RuntimeError(message)

            # next, the parameters have to be valid BAYESN params
            if paramname not in ALLOWED_BAYESN_PARAMS:
                message = f'Do not recognize key {key} in ARGLIST.\n'
                message += f'Keys must specify a BayeSN light curve parameter f{ALLOWED_BAYESN_PARAMS}.'
                raise RuntimeError(message)

            # we need to keep a track of which distribution parameters are set
            temp = self.bayesn_param_config.get(paramname, None)
            if temp is None:
                temp = []
            temp.append(keysuffix)
            self.bayesn_param_config[paramname] = temp

            # finally you either get a single float or a sequence of floats
            # to specify the distributions for each parameter
            try:
                temp = float(val)
                if paramname == 'RANGE':
                    message = f'Key {key} specified but only found a single value.'
                    raise RuntimeError(message)
            except ValueError:
                temp = [float(v) for v in val.split()]
                if keysuffix != 'RANGE':
                    message = f'Key {key} specified and has multiple values, but only accepts a single value.'
                    raise RuntimeError(message)
                if len(temp) !=2 :
                    message = f'Range key {key} specified but can only have two values.'
                    raise RuntimeError(message)
                if temp[0] >= temp[1]:
                    message = f'Range key {key} specified but lower limit is >= uppler limit.'
                    raise RuntimeError(message)
            self.bayesn_distrib_params[key] = temp
        print('Parsed the following BayeSN distribution parameters')
        print(self.bayesn_distrib_params)


    def setup_bayesn_param_distributions(self):
        """
        Parse the BayeSN parameter distributions specified in
        bayesn.params or as a string in the ARGLIST in the SNANA .INPUT
        file and configure distribution objects

        These distribution objects are sampled to draw parameters for
        generating the light curves

        Currently the options are
        set a scalar (PEAK or PEAK and RANGE),
        set a Uniform distribution  (RANGE, but no PEAK or SIGMA)
        set a Gaussian (PEAK AND SIGMA)
        or set a truncated Gaussian (PEAK, SIGMA and RANGE)

        Distributions are saved as objects in the dictionary
        bayesn_distribs as
        <PARAM>: <function object that can be called>
        """
        self.bayesn_distribs = {}
        for param, param_distrib_names in self.bayesn_param_config.items():
            if len(param_distrib_names) == 1:
                # if only one parameter is specified then
                # it had better be peak or range
                    distrib_name = param_distrib_names[0]
                    key = f'GEN{distrib_name}_{param}'
                    val = self.bayesn_distrib_params[key]
                    if distrib_name == 'RANGE':
                        # We have a range, so set a random uniform distribution
                        llim, ulim = val
                        scale = ulim - llim
                        distrib = st.uniform(loc=llim, scale=scale)
                        self.bayesn_distribs[param] = distrib

                    elif distrib_name == 'PEAK':
                        # we have only a peak - just set this as a scalar
                        self.setParVals_BAYESN(param=float(val))
            else:
                # if we have more than one parameter specified
                # then we are specifying a Gaussian of some sort
                peak  = self.bayesn_distrib_params.get(f'GENPEAK_{param}', None)
                sigma = self.bayesn_distrib_params.get(f'GENSIGMA_{param}', None)
                llim, ulim  = self.bayesn_distrib_params.get(f'GENRANGE_{param}', (None, None))

                # if it's PEAK AND RANGE - OK IF PEAK WITHIN RANGE
                if llim is not None and ulim is not None:
                    if peak is not None:
                        if llim >= peak and peak >= ulim:
                            message = f'{param} range {llim}--{ulim} peak {peak} specified but peak is outside range'
                            raise RuntimeError(message)

                # if it's SIGMA AND RANGE - NOT OK
                if peak is None and sigma is not None:
                        message = f'{param} peak {peak} not specified but sigma {sigma} is specified - invalid'
                        raise RuntimeError(message)

                # if peak is specified and sigma is not, we've either got a range and checked
                # or we've not got a range and whatever - it's still a valid fixed value
                if sigma is None and peak is not None:
                    self.setParVals_BAYESN(param=float(peak))

                if peak is not None and sigma is not None:
                    if llim is not None and ulim is not None:
                        # if it's all three - OK IF PEAK WITHIN RANGE
                        # we've check that already
                        a, b =  (llim - peak) / sigma, (ulim - peak) / sigma
                        distrib = st.truncnorm(a, b, loc=peak, scale=sigma)
                    else:
                        # if it's PEAK AND SIGMA - OK
                        distrib = st.norm(loc=peak, scale=sigma)
                    self.bayesn_distribs[param] = distrib
            # end logic for number of parameters set
        # end loop over parameters


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
            useful_pars = {x:hostpars[i] for i, x in enumerate(self.host_param_names)}
            for param, distrib in self.bayesn_distribs.items():
                useful_pars[param] = distrib.rvs()

            eta = np.random.normal(0, 1, len(self._bayesn_components["L_Sigma_epsilon"]))
            epsilon_vec = np.dot(self._bayesn_components["L_Sigma_epsilon"], eta)
            for i in range(self._nepsilon):
                useful_pars[f'EPSILON{i:02d}'] = epsilon_vec[i]
            if self.verbose:
                print(useful_pars, 'set')

            self.setParVals_BAYESN(**useful_pars)


        #ST: Three lines originally from GSN
        ind =  np.abs(self.phase - trest).argmin()
        ind_flux = ind*self.wavelen
        flux = self.flux[ind_flux:ind_flux+self.wavelen]

        ########## NEW CODE FROM ST BELOW HERE (FOR GUIDANCE) ##########


        #ST: Computes matrices that do interpolation
        #    Probably can't be precomputed
        #    Assumes that `trest` is a float, and `self.wave` is a 1D
        #    list or numpy array of rest frame wavelengths
        J_t =  spline_coeffs_irr([trest], self._bayesn_components["tau_knots"], self.KD_t).T

        #ST: Computes host extinction
        #    This assumes we can use the Kyle Barbary extinction.py package
        #    If we can't, I have the necessary code for this
        #    Note: This may need `self.wave` to be converted to a numpy
        #    array, if it is a list.
        R_host = extinction.fitzpatrick99(self.wave, self.parameter_values["AV"], self.parameter_values["RV"])


        #ST: Computes the spline knots
        epsilon_matrix = np.zeros(self._bayesn_components["W0"].shape)
        epsilon_vector = [self.parameter_values[f'EPSILON{i:02d}'] for i in range(self._nepsilon)]
        epsilon_matrix[1:-1,:] = np.reshape(epsilon_vector, epsilon_matrix[1:-1,:].shape, order="F")
        W = self._bayesn_components["W0"] + self.parameter_values["THETA1"]*self._bayesn_components["W1"] + epsilon_matrix

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


def invKD_irr(x):
	"""
	Compute K^{-1}D for a set of spline knots.

	For knots y at locations x, the vector, y'' of non-zero second
	derivatives is constructed from y'' = K^{-1}Dy, where K^{-1}D
	is independent of y, meaning it can be precomputed and reused for
	arbitrary y to compute the second derivatives of y.

	Parameters
	----------
	x : :py:class:`numpy.array`
		Numpy array containing the locations of the cubic spline knots.

	Returns
	-------
	KD : :py:class:`numpy.array`
		y independednt matrix whose product can be taken with y to
		obtain a vector of second derivatives of y.
	"""
	n = len(x)

	K = np.zeros((n-2,n-2))
	D = np.zeros((n-2,n))

	K[0,0:2] = [(x[2] - x[0])/3, (x[2] - x[1])/6]
	K[-1, -2:n-2] = [(x[n-2] - x[n-3])/6, (x[n-1] - x[n-3])/3]

    # should be able to vectorize this - GN 05/16/22
	for j in np.arange(2,n-2):
		row = j - 1
		K[row, row-1:row+2] = [(x[j] - x[j-1])/6, (x[j+1] - x[j-1])/3, (x[j+1] - x[j])/6]
	for j in np.arange(1,n-1):
		row = j - 1
		D[row, row:row+3] = [1./(x[j] - x[j-1]), -(1./(x[j+1] - x[j]) + 1./(x[j] - x[j-1])), 1./(x[j+1] - x[j])]

	M = np.zeros((n,n))
	M[1:-1, :] = np.linalg.solve(K,D)
	return M


def cartesian_prod(x, y):
	"""
	Compute cartesian product of two vectors.

	Parameters
	----------
	x : :py:class:`numpy.array`
		First vector.
	x : :py:class:`numpy.array`
		Second vector.

	Returns
	-------
	z : :py:class:`numpy.array`
		Cartesian product of x and y.
	"""
	n_x = len(x)
	n_y = len(y)
	return np.array([np.repeat(x,n_y),np.tile(y,n_x)]).T


def spline_coeffs_irr(x_int, x, invkd, allow_extrap=True):
	"""
	Compute a matrix of spline coefficients.

	Given a set of knots at x, with values y, compute a matrix, J, which
	can be multiplied into y to evaluate the cubic spline at points
	x_int.

	Parameters
	----------
	x_int : :py:class:`numpy.array`
		Numpy array containing the locations which the output matrix will
		interpolate the spline to.
	x : :py:class:`numpy.array`
		Numpy array containing the locations of the spline knots.
	invkd : :py:class:`numpy.array`
		Precomputed matrix for generating second derivatives. Can be obtained
		from the output of ``invKD_irr``.
	allow_extrap : bool
		Flag permitting extrapolation. If True, the returned matrix will be
		configured to extrapolate linearly beyond the outer knots. If False,
		values which fall out of bounds will raise ValueError.

	Returns
	-------
	J : :py:class:`numpy.array`
		y independednt matrix whose product can be taken with y to evaluate
		the spline at x_int.
	"""
	n_x_int = len(x_int)
	n_x = len(x)
	X = np.zeros((n_x_int,n_x))

	if not allow_extrap and ((max(x_int) > max(x)) or (min(x_int) < min(x))):
		raise ValueError("Interpolation point out of bounds! " +
			"Ensure all points are within bounds, or set allow_extrap=True.")

	for i in range(n_x_int):
		x_now = x_int[i]
		if x_now > max(x):
			h = x[-1] - x[-2]
			a = (x[-1] - x_now)/h
			b = 1 - a
			f = (x_now - x[-1])*h/6.0

			X[i,-2] = a
			X[i,-1] = b
			X[i,:] = X[i,:] + f*invkd[-2,:]
		elif x_now < min(x):
			h = x[1] - x[0]
			b = (x_now - x[0])/h
			a = 1 - b
			f = (x_now - x[0])*h/6.0

			X[i,0] = a
			X[i,1] = b
			X[i,:] = X[i,:] - f*invkd[1,:]
		else:
			q = np.where(x[0:-1] <= x_now)[0][-1]
			h = x[q+1] - x[q]
			a = (x[q+1] - x_now)/h
			b = 1 - a
			c = ((a**3 - a)/6)*h**2
			d = ((b**3 - b)/6)*h**2

			X[i,q] = a
			X[i,q+1] = b
			X[i,:] = X[i,:] + c*invkd[q,:] + d*invkd[q+1,:]

	return X


def main():
    mySED=gensed_BAYESN('$SNDATA_ROOT/models/bayesn/BAYESN.M20',2,[],'z,AGE,ZCMB,METALLICITY')

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(111)
    mySED.setParVals_BAYESN(THETA1=0., AV=0.1 , RV=3.1 , DELTAM=0.  , TMAX=0.)
    flux = mySED.fetchSED_BAYESN(0, new_event=2)
    ax1.plot(mySED.wave, np.array(flux), 'g-')

    mySED.setParVals_BAYESN(THETA1=2., AV=0.3 , RV=3.1 , DELTAM=0.  , TMAX=0.)
    flux = mySED.fetchSED_BAYESN(0, new_event=2)
    ax1.plot(mySED.wave, np.array(flux), 'r-')

    mySED.setParVals_BAYESN(THETA1=-1.7, AV=0. , RV=3.1 , DELTAM=0.5 , TMAX=0.)
    flux = mySED.fetchSED_BAYESN(0, new_event=2)
    ax1.plot(mySED.wave, np.array(flux), 'b-')

    ax1.plot(mySED.wave, np.array(flux))
    fig.savefig('test.png')

    fig2 = plt.figure(figsize=(10, 20))
    ax2 = fig2.add_subplot(111)
    ind = (mySED.wave >= 3000) & (mySED.wave <= 7000)
    for trest in np.linspace(-10, 40, 20):

        mySED.setParVals_BAYESN(THETA1=-1.7, AV=0. , RV=3.1 , DELTAM=0. ,  TMAX=0.)
        flux = mySED.fetchSED_BAYESN(trest, new_event=2)
        ax2.plot(mySED.wave[ind], np.array(flux)[ind]/np.array(flux)[ind].max() + trest, 'b-')

        mySED.setParVals_BAYESN(THETA1=0., AV=0.1 , RV=3.1 , DELTAM=0. , TMAX=0.)
        flux = mySED.fetchSED_BAYESN(trest, new_event=2)
        ax2.plot(mySED.wave[ind], np.array(flux)[ind]/np.array(flux)[ind].max() + trest, 'g-')

        mySED.setParVals_BAYESN(THETA1=2., AV=0.3 , RV=3.1 , DELTAM=0. , TMAX=0.)
        flux = mySED.fetchSED_BAYESN(trest, new_event=2)
        ax2.plot(mySED.wave[ind], np.array(flux)[ind]/np.array(flux)[ind].max() + trest, 'r-')
    fig2.savefig('phase_sequence.png')


if __name__=='__main__':
    main()
