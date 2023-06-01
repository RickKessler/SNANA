"""
This is the model file for AGN object. We adopted damped random walk model.

Description
The model needs time series provided by SNANA. The black hole mass, eddington ratio, Mi, wavelength range, rng(random seed) are provided by the model.
The model will return flux at provided wavelength at each time moment.
For debug purpose, you can run this model as python script without SNANA. See main() function for details.

SNANA Inputs:
cl_prob: probability that changing look behavior will happen within 1 year. Currently, changing look behavior happens
        at most once in one AGN event.

SNANA Outputs:
Mi: characteristic i-band absolute magnitude K-corrected to z=2
M_BH: Black hole mass in unit of solar mass
edd_ratio: Eddington ratio
edd_ratio2: Eddington ratio after t_transition for CLAGN. If it's regular AGN, edd_ratio2 is the same as edd_ratio.
t_transition: rest frame time of changling look transition. We output this in both regular AGN and CLAGN cases.
            But no transition takes place in regular AGN.
            In oberver's frame, trasition date is PEAKMJD + t_transition * (1 + z)

To get these outputs, you can simply add the following to .INPUT file:
SIMGEN_DUMP:  17
CID GENTYPE  SNTYPE  NON1A_INDEX  GENZ
LIBID RA DECL MWEBV MU PEAKMJD
M_BH Mi edd_ratio edd_ratio2 t_transition cl_flag

Input parameters are given by AGN.INFO, located at path given by "GENMODEL AGN PATH". Example of AGN.INFO file:
CONFIG:
  cl_prob: 1e-2
You can overwrite parameters from AGN.INFO in input file using ARGLIST statement. 

References:
Mi derivation for AGN:
    Shen et al., 2013
    https://adsabs.harvard.edu/full/2013BASI...41...61S
Standard Disk Model:
    Lipunova, G., Malanchev, K., Shakura, N. (2018). Page 33 for the main equation
    DOI https://doi.org/10.1007/978-3-319-93009-1_1
Time scale (tau) and SF_inf parameters:
    Suberlak et al. 2021
    DOI 10.3847/1538-4357/abc698
Eddington ratio distribution:
        https://github.com/burke86/imbh_forecast/blob/master/var.ipynb
        https://ui.adsabs.harvard.edu/abs/2019ApJ...883..139S/abstract
        parameters from this paper:
        https://iopscience.iop.org/article/10.3847/1538-4357/aa803b/pdf
"""
from pathlib import Path
import numpy as np
from astropy import constants, units
from scipy import integrate
from scipy.interpolate import UnivariateSpline
import re
from gensed_base import gensed_base
import yaml


M_sun = constants.M_sun.cgs.value
c = constants.c.cgs.value
pc = (1 * units.pc).to_value(units.cm)
G = constants.G.cgs.value
sigma_sb = constants.sigma_sb.cgs.value
h = constants.h.cgs.value
k_B = constants.k_B.cgs.value
year = units.year.to(units.day)


def search_nearest(a, x):
    """Find index of element of `a` closest to x

    Parameters
    ----------
    a : np.array
        Sorted array
    x : number
        Value to loop up

    Returns
    -------
    int
        Valid index of a
    """
    idx = np.searchsorted(a, x)
    if idx == 0:       # x < a[0]
        return 0
    if idx == a.size:  # x >= a[-1]
        return a.size - 1
    if x - a[idx - 1] < a[idx] - x:
        return idx - 1
    return idx


class DistributionSampler:
    """
    A class to build a sampler by inverse transformation sampling

    ...

    Attributes
    ----------
    None

    Methods
    -------
    inv_cdf(x, y)
        Returns inverse cdf function
    inv_trans_sampling(self, sampling_size=1)
        Returns random samples by inverse transformation sampling
    """

    def __init__(self, x, y, rng=None):
        """
        Parameters
        ----------
        x : ndarray[float64]
            The random variable
        y : str
            The probability density function of x
        rng : random generator object, optional
            The initialized generator object for inverse transform sampling
        """

        self.rng = np.random.default_rng(rng)
        self.inv_cdf_spline = self.inv_cdf(x, y)

    @staticmethod
    def inv_cdf(x, y):
        """
        Return inverse cumulative distribution function

        Parameters
        ----------
        x : ndarray[float64]
            The random variable
        y : str
            The probability density function of x
        """

        pdf_spline = UnivariateSpline(x=x, y=y, s=0, k=3)
        cdf_spline = pdf_spline.antiderivative()
        assert cdf_spline(x[0]) == 0
        norm = cdf_spline(x[-1])
        inv_cdf_spline = UnivariateSpline(x=cdf_spline(x) / norm, y=x, s=0, k=3)

        return inv_cdf_spline

    def inv_trans_sampling(self, sampling_size=1):
        """
        Return samples from inverse transformation sampling

        Parameters
        ----------
        sampling_size : int
            Desired sampling size
        """

        r = self.rng.random(int(sampling_size))
        return self.inv_cdf_spline(r)


class AGN:
    """
    A class for AGN object

    ...

    Arguments
    ---------
    t0: float
        initial time moment [day]
    Mi: float
        characteristic i-band absolute magnitude K-corrected to z=2
    M_BH: float
        black hole mass, in unit of g
    lam: np.ndarray
        wavelength
    edd_ratio: float
        Eddington ratio
    rng:
        random seed

    Attributes
    ----------
    lam: ndarray
        array of rest frame wavelength
    t0: float
        initial time moment
    rng: random seed generator object
        random seed
    ME_dot: float
        Accretion rate at Eddington luminosity
    MBH_dot: float
        Black hole accretion rate
    Fnu_average: ndarray[float64]
        Baseline for flux density Fnu
    tau: ndarray
        Timescale for damped random walk model at each wavelength
    sf_inf:  ndarray
        Structure Function (infinity) at each wavelength
    t: float
        current time moment for damped random walk model
    """
    def __init__(self, t0: float, M_BH: float, lam: np.ndarray, edd_ratio: float, rng):
        self.lam = np.asarray(lam)
        self.t0 = t0
        self.rng = np.random.default_rng(rng)

        self.ME_dot = self.find_ME_dot(M_BH)
        self.MBH_dot = self.find_MBH_dot(self.ME_dot, edd_ratio)
        self.Fnu_average = 2 * self.find_Fnu_average_standard_disk(self.MBH_dot, self.lam, M_BH) # quick fix to double the baseline Fnu

        L_bol = self.find_L_bol(edd_ratio, M_BH)
        self.Mi = self.find_Mi(L_bol)
        self.tau = self.find_tau_v(self.lam, self.Mi, M_BH)
        self.sf_inf = self.find_sf_inf(self.lam, self.Mi, M_BH)

        self.t = t0
        self.delta_m = self._random() * self.sf_inf

    def step(self, t):
        """
        time step of damped random walk model. t must be larger than self.t
        """
        dt = t - self.t
        self.t = t

        self.delta_m = (
                self.delta_m * np.exp(-dt / self.tau)
                + self.sf_inf * np.sqrt(1 - np.exp(-2 * dt / self.tau)) * self._random()
        )

    def __call__(self, t):
        self.step(t)
        return self.Fnu

    @staticmethod
    def find_L_bol(edd_ratio, M_BH):
        """
        Input: M_BH in unit of g.
        Return L_bol in erg/s
        """
        return edd_ratio * 1.26e38 * M_BH /M_sun

    @staticmethod
    def find_Mi(L_bol):
        """
        Input: L_bol in erg/s
        Outbout: Mi
        """
        # from Shen et al., 2013
        # https://adsabs.harvard.edu/full/2013BASI...41...61S
        return 90 - 2.5 * np.log10(L_bol)

    @property
    def Fnu(self):
        """
        Correct Fnu with its baseline (average value)
        Output: Fnu: spectral density along the frequency axes
        """
        return 10 ** (-0.4 * self.delta_m) * self.Fnu_average

    def _random(self):
        return self.rng.normal(size=self.lam.size)

    @staticmethod
    def find_ME_dot(M_BH):
        """
        Input: BH mass in g
        return: Accretion rate at Eddington luminosity in g/s
        """
        return 1.4e18 * M_BH / M_sun

    @staticmethod
    def find_MBH_dot(ME_dot, eddington_ratio):
        """
        Input:
        ME_dot: Accretion rate at Eddington luminosity in g/s
        eddington ratio
        Return:
        MBH_dot: Black hole accretion rate
        """
        return ME_dot * eddington_ratio

    @staticmethod
    def T_0(M, Mdot, r_in):
        """
        calculate T0 based on standard disk model
        input:
        M: mass of the gravitating centre
        Mdot: accretion rate at previous time step
        r_in: the inner radius of the accretion disc
        output:
        T_0: Effective temperature at r0. Same as the maximum effective temperature at the disc surface (Tmax)
        """
        # Lipunova, G., Malanchev, K., Shakura, N. (2018). Page 33 for the main equation
        # DOI https://doi.org/10.1007/978-3-319-93009-1_1

        return (2 ** (3 / 4) * (3 / 7) ** (7 / 4) * (G * M * Mdot / (np.pi * sigma_sb * r_in ** 3)) ** (1 / 4))

    @staticmethod
    def r_0(r_in):
        """
        return r0 based on standard disk model

        input: r_in: the inner radius of the accretion disc
        output: r0:  the initial radius of the ring
        """
        # Lipunova, G., Malanchev, K., Shakura, N. (2018). Page 33 for the main equation
        # DOI https://doi.org/10.1007/978-3-319-93009-1_1

        return ((7 / 6) ** 2 * r_in)

    @staticmethod
    def x_fun(nu, T0, r, r0):
        """
        calculate variable of integration x

        input:
        nu: frequency
        T0: Effective temperature at r0. Same as the maximum effective temperature at the disc surface (Tmax)
        r: radius of the accretion disc
        r0: the initial radius of the ring

        return: variable of integration x
        """
        # Lipunova, G., Malanchev, K., Shakura, N. (2018). Page 33 for the main equation
        # DOI https://doi.org/10.1007/978-3-319-93009-1_1
        return (h * nu / (k_B * T0) * (r / r0) ** (3 / 4))

    def find_flux_standard_disk(self, Mdot, nu, rin, i, d, M):
        """
        function to calculate flux based on standard disk model

        input:
        Mdot: accretion rate at previous time step
        nu: frequency
        rin: the inner radius of the accretion disc
        i: inclination
        d: distance
        M: mass of the gravitating centre
        output:
        flux at given time step
        """
        # Lipunova, G., Malanchev, K., Shakura, N. (2018). Page 33 for the main equation
        # DOI https://doi.org/10.1007/978-3-319-93009-1_1
        T0 = self.T_0(M, Mdot, rin)
        r0 = self.r_0(rin)
        # large x in exponetial causes overflow, but 1/inf is zero.
        with np.errstate(over='ignore'):
            fun_integr = lambda x: (x ** (5 / 3)) / np.expm1(x)
            integ, inte_err = integrate.quad(fun_integr, 1e-6, np.inf)

        return ((16 * np.pi) / (3 * d ** 2) * np.cos(i) * (k_B * T0 / h) ** (8 / 3) * h * (nu ** (1 / 3)) / (c ** 2) * (
                r0 ** 2) * integ)

    def find_Fnu_average_standard_disk(self, MBH_dot, lam, M_BH):
        flux_av = self.find_flux_standard_disk(MBH_dot, c / lam, rin=1, i=0, d=10 * pc,
                                               M=M_BH)  # SNANA required flux observed from 10 pc
        return flux_av

    @staticmethod
    def find_tau_v(lam, Mi=-23, M_BH=1e9 * M_sun):
        """
        Return timescale for DRW model.
        Input frequency v in Hz, i band magnitude (default is -23), Black hole mass in g (defalt is 10^9 solar mass).
        Return timescale in s.
        """
        # Equation and parameters for A, B, C, D adopted from Suberlak et al. 2021
        # DOI 10.3847/1538-4357/abc698

        A = 2.4
        B = 0.17
        C = 0.03
        D = 0.21
        return 10 ** (A + B * np.log10(lam / (4000e-8))
                      + C * (Mi + 23) + D * np.log10(M_BH / (1e9 * M_sun)))  # e-8 angstrom

    @staticmethod
    def find_sf_inf(lam, Mi=-23, M_BH=1e9 * M_sun):
        """
        Return Structure Function at infinity in mag.
        Input frequency in Hz, i band magnitude (default is -23), Black hole mass in g (defalt is 10^9 solar mass).
        """
        # Equation and parameters for A, B, C, D adopted from Suberlak et al. 2021
        # DOI 10.3847/1538-4357/abc698
        A = -0.51
        B = -0.479
        C = 0.13
        D = 0.18

        return 10 ** (A + B * np.log10(lam / (4000e-8))
                      + C * (Mi + 23) + D * np.log10(M_BH / (1e9 * M_sun)))


class gensed_AGN(gensed_base):
    def __init__(self, PATH_VERSION, OPTMASK, ARGLIST, HOST_PARAM_NAMES):
        print('__init__', flush=True)
        print(HOST_PARAM_NAMES)
        print(PATH_VERSION)
        self.agn1 = None
        self.agn2 = None
        self.trest = None
        self.t_transition = None
        self.sed = None
        self.sed_Fnu = None

        self.parse_param(arglist=ARGLIST, path_version = PATH_VERSION)
        self.cl_prob_event = None
        self.cl_flag = False

        self.rng = np.random.default_rng(0)

        self.wavelen = 100
        self.wave = np.logspace(np.log10(100e-8), np.log10(20000e-8), self.wavelen)

        self.edd_ratio = None
        self.edd_ratio2 = None

        self.log_lambda_min = -8.5
        self.log_lambda_max = 0.5
        self.nbins = 1000

        self.M_BH = None
        self.Mi = None

        lambda_ = np.logspace(self.log_lambda_min, self.log_lambda_max, self.nbins + 1)
        xi_blue = self.ERDF(lambda_Edd=lambda_, rng=self.rng)
        self.ERDF_spline = DistributionSampler(lambda_, xi_blue, rng=self.rng)
        print('gensed_AGN model initialized')

    def parse_param(self, arglist, path_version):
        path_version = Path(path_version)
        if not path_version.exists():
            raise ValueError(f'genmodel AGN path does not exist: {path_version}')
        with open(path_version/'AGN.INFO') as file:
            info = yaml.load(file, yaml.Loader)
        config = info['CONFIG']

        args = arglist.split()
        if len(args) % 2 != 0:
            raise ValueError(f'Invalid argument list, please provide key-value pairs instead: "ARGLIST:{arglist}"')
        args = dict(zip(args[::2], args[1::2]))

        config.update(args)

        self.ranseed = int(config.pop('RANSEED', 0))
        self.cl_prob = float(config.pop('cl_prob', 0.0))
        if len(config) != 0:
            raise ValueError(f'ARGLIST or {path_version/"AGN.INFO"} has unknown arguments: {" ".join(config)}')


    def Fnu_to_Flamb(self, Fnu):
        """
        convert from Fnu [erg/s/cm^2/Hz] to Flamb[erg/s/cm^2/AA]
        """
        return Fnu * c / self.wave ** 2 * 1e-8

    def Flamb_to_Fnu(self, Flamb):
        """
        convert from to Flamb[erg/s/cm^2/AA] to Fnu [erg/s/cm^2/Hz]
        """
        return Flamb * self.wave ** 2 * 1e8 / c


    @staticmethod
    def ERDF(lambda_Edd, galaxy_type='Blue', rng=None):
        """
        Eddington Ratio Distribution Function for blue galaxies (radiatively-efficient, less massive)
        """
        # https://github.com/burke86/imbh_forecast/blob/master/var.ipynb
        #
        # https://ui.adsabs.harvard.edu/abs/2019ApJ...883..139S/abstract
        # parameters from this paper:
        # https://iopscience.iop.org/article/10.3847/1538-4357/aa803b/pdf
        rng = np.random.default_rng(rng)
        # Lbr = 10**38.1 lambda_br M_BH_br
        # 10^41.67 = 10^38.1 * 10^x * 10^10.66
        if galaxy_type == 'Red':
            xi = 10 ** -2.13
            lambda_br = 10 ** rng.normal(-2.81, np.mean([0.22, 0.14]))
            delta1 = rng.normal(0.41 - 0.7, np.mean([0.02, 0.02]))  # > -0.45 won't affect LF
            delta2 = rng.normal(1.22, np.mean([0.19, 0.13]))

        if galaxy_type == 'Blue':
            xi = 10 ** -1.65
            lambda_br = 10 ** rng.normal(-1.84, np.mean([0.30, 0.37]))
            delta1 = rng.normal(0.471 - 0.7, np.mean([0.02, 0.02]))  # > -0.45 won't affect LF
            delta2 = rng.normal(2.53, np.mean([0.68, 0.38]))

        return xi * ((lambda_Edd / lambda_br) ** delta1 + (lambda_Edd / lambda_br) ** delta2) ** -1

    @staticmethod
    def M_BH_sample(rng):
        logMH_min = 7
        logMBH_max = 9
        return 10 ** (rng.uniform(logMH_min, logMBH_max))


    def fetchSED_NLAM(self):
        """
        Returns the length of the wavelength vector
        """
        print('fetchSED_NLAM', flush=True)
        return self.wavelen


    def prepEvent(self, trest, external_id, hostparams):
        """
        generate the full event, saved to self.sed, and saved to self.trest
        """
        self.edd_ratio = self.ERDF_spline.inv_trans_sampling(sampling_size=1)
        self.M_BH = self.M_BH_sample(self.rng)  # M_BH in unit of M_sun

        self.trest = np.unique(trest)
        duration = self.trest[-1] - self.trest[0]
        self.cl_prob_event = 1 - (1-self.cl_prob)**(duration/year)
        # select CL transition time moment. We use base AGN model for trest[0] to t_transition, and we use model
        # with higher eddington ratio for t_transition to trest[-1]
        self.t_transition = self.rng.uniform(trest[0], trest[-1])
        trans_idx = np.searchsorted(self.trest, self.t_transition)
        trest1 = self.trest[:trans_idx]
        trest2 = self.trest[trans_idx:]

        self.agn1 = AGN(t0=trest1[0], M_BH=self.M_BH * M_sun, lam=self.wave, edd_ratio=self.edd_ratio,
                       rng=self.rng)
        self.Mi = self.agn1.Mi
        if self.rng.random() < self.cl_prob_event:
            self.edd_ratio2 = self.rng.uniform(5, 30) * self.edd_ratio
            self.cl_flag = True
            # print('CLAGN')
        else:
            self.edd_ratio2 = self.edd_ratio
            self.cl_flag = False
            # print('regular AGN')

        self.agn2 = AGN(t0=self.t_transition, M_BH=self.M_BH * M_sun, lam=self.wave,
                        edd_ratio=self.edd_ratio2, rng=self.rng)

        # print(self.agn2.tau[0]/self.agn1.tau[0])
        # print(self.agn2.sf_inf[0] / self.agn1.sf_inf[0])

        # initial SED is randomly sampled.
        sed_lam = [self.Fnu_to_Flamb(self.agn1.Fnu)]
        # Do DRW stepping for AGN1
        for t in trest1[1:]:
            self.agn1.step(t)
            sed_lam.append(self.Fnu_to_Flamb(self.agn1.Fnu))
        # Do one more stepping for AGN1 for smooth transition to AGN2. We do not append this to F_lam intentionally.
        self.agn1.step(self.t_transition)
        # In AGN2 constructor, we initialize delta_m randomly.
        # Here instead, we set it to produce the same Fnu as the last step of AGN1
        self.agn2.delta_m = self.agn1.delta_m + 2.5 * np.log10(self.agn2.Fnu_average / self.agn1.Fnu_average)
        # Do DRW stepping for ANG2
        for t in trest2:
            self.agn2.step(t)
            sed_lam.append(self.Fnu_to_Flamb(self.agn2.Fnu))
        self.sed = np.stack(sed_lam)

    def test_AGN_flux(self, trest):
        """
        for debug purpose. Only used in main() function.
        """
        self.trest = np.unique(trest)
        duration = self.trest[-1] - self.trest[0]
        self.cl_prob_event = 1 - (1-self.cl_prob)**(duration/year)
        # select CL transition time moment. We use base AGN model for trest[0] to t_transition, and we use model
        # with higher eddington ratio for t_transition to trest[-1]
        self.t_transition = self.rng.uniform(trest[0], trest[-1])
        trans_idx = np.searchsorted(self.trest, self.t_transition)
        trest1 = trest[:trans_idx]
        trest2 = trest[trans_idx:]
        print(self.t_transition)
        print(trest[trans_idx])

        self.agn1 = AGN(t0=trest1[0], M_BH=self.M_BH * M_sun, lam=self.wave, edd_ratio=self.edd_ratio,
                       rng=self.rng)
        self.Mi = self.agn1.Mi
        if self.rng.random() < self.cl_prob_event:
            self.edd_ratio2 = self.rng.uniform(5, 30) * self.edd_ratio
            self.cl_flag = True
            print('CLAGN')
        else:
            self.edd_ratio2 = self.edd_ratio
            self.cl_flag = False
            print('regular AGN')

        self.agn2 = AGN(t0=self.t_transition, M_BH=self.M_BH * M_sun, lam=self.wave,
                        edd_ratio=self.edd_ratio2, rng=self.rng)

        print(self.agn2.tau[0]/self.agn1.tau[0])
        print(self.agn2.sf_inf[0] / self.agn1.sf_inf[0])

        # initial SED is randomly sampled.
        sed_lam = [self.Fnu_to_Flamb(self.agn1.Fnu)]
        # Do DRW stepping for AGN1
        for t in trest1[1:]:
            self.agn1.step(t)
            sed_lam.append(self.Fnu_to_Flamb(self.agn1.Fnu))
        # Do one more stepping for AGN1 for smooth transition to AGN2. We do not append this to F_lam intentionally.
        self.agn1.step(self.t_transition)
        # In AGN2 constructor, we initialize delta_m randomly.
        # Here instead, we set it to produce the same Fnu as the last step of AGN1
        self.agn2.delta_m = self.agn1.delta_m + 2.5 * np.log10(self.agn2.Fnu_average / self.agn1.Fnu_average)
        # Do DRW stepping for ANG2
        for t in trest2:
            self.agn2.step(t)
            print(self.agn2.delta_m[61])
            sed_lam.append(self.Fnu_to_Flamb(self.agn2.Fnu))
        self.sed = np.stack(sed_lam)


    def fetchSED_LAM(self):
        """
        Returns the wavelength vector
        """
        wave_aa = self.wave * 1e8
        return wave_aa

    def fetchSED(self, trest, maxlam, external_id, new_event, hostparams):
        # We get observation created from prep_event()

        idx = search_nearest(self.trest, trest)
        # Check if given time value is close enough to what we've found
        # 1e-4 days is arbitrary small value
        if abs(self.trest[idx] - trest) > 1e-4:
            raise ValueError(f"trest = {trest:.6f} was not providen in prepEvent()")
        return self.sed[idx]

    def fetchParNames(self):
        return ['M_BH', 'Mi', 'edd_ratio', 'edd_ratio2', 't_transition', 'cl_flag']

    def fetchParVals(self, varname):
        return getattr(self, varname)

def main():
    """
    We provide 3 sets of input values to generate ligthcurve for debugging purpose. 2 of them are realistic objects with dbID available
    Code to plot eddington ratio distribution, spectrum available at the bottom of main(), comment out to use them
    """
    import matplotlib.pyplot as plt
    import astropy
    from astropy.cosmology import Planck18

    mySED = gensed_AGN('$SNDATA_ROOT/models/bayesn/BAYESN.M20',2,'RANSEED 12345 cl_prob 0.99','z,AGE,ZCMB,METALLICITY')
    trest = np.arange(1, 1000, 0.1)

    # test set 1  for obj dbID=2534406
    # z = 0.805899978
    L_bol = 10**(45.21136675358238)
    mySED.M_BH = 10 ** (8.564795254)
    mySED.edd_ratio = L_bol/ (1.26e38 * mySED.M_BH)
    mySED.rng = np.random.default_rng(0)
    # mySED.Mi = mySED.find_Mi(L_bol)

    ##test set 2 for obj dbID=8442
    ## z = 2.43
    # L_bol = 10**(46.61)
    # mySED.M_BH = 10 ** 9.09
    # mySED.edd_ratio = L_bol/ (1.26e38 * mySED.M_BH)
    # mySED.rng = np.random.default_rng(0)
    # mySED.Mi = mySED.find_Mi(L_bol)

    # #test set 3
    # mySED.Mi = -23
    # mySED.M_BH = 1e9 # in unit of M_sun
    # mySED.rng = np.random.default_rng(0)
    # mySED.edd_ratio = 0.1


    mySED.test_AGN_flux(trest)
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(111)

    flux_firstWave = []
    sed_list = - 2.5 * np.log10(mySED.Flamb_to_Fnu(mySED.sed)) - 48.5
    # sed_list = - 2.5 * np.log10(list(mySED.sed_Fnu.values())) - 48.5
    sed_list = sed_list + astropy.coordinates.Distance(z=0.805899978).distmod.value  # adjust z for test set 1 and 2

    flux_firstWave = sed_list[:, 61]  #wave[61] closest to 4770/(1+ 0.805899978) armstrong(corrected SDSS g band)  dbID=2534406, comment out for test set 1
    #flux_firstWave.append(sed_list[:][54])  #wave[54] closest to 6156/(1+2.43) armstrong(corrected ps1 r band), comment out for test set 2
    #flux_firstWave.append(sed_list[:][82])  #wave[82] = 8000 armstrong(i-band), comment out for test set 3
    ax1.vlines(mySED.t_transition*(1+0.805899978), 19, 21)
    print(mySED.agn1.tau[61])

    ax1.plot(trest*(1+0.805899978), flux_firstWave, 'g-')
    ax1.invert_yaxis()
    plt.show()
    fig.savefig('debug_apparentmag_2534406.png')


    # ## test for ERDF range
    # mySED.edd_ratio = mySED.ERDF_spline.inv_trans_sampling(sampling_size=1000)
    # plt.hist(mySED.edd_ratio, bins=100)
    # plt.yscale('log')
    # plt.title('Eddinton ratio distribution')
    # plt.show()


    ### plot spectrum
    # F_spectrum = mySED.sed_Fnu[10]
    # F_spectrum2 = mySED.sed_Fnu[20]
    #
    # ax1.plot(mySED.wave, F_spectrum)
    # ax1.plot(mySED.wave, F_spectrum2)

if __name__ == '__main__':
    main()
