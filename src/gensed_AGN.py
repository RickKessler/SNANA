import numpy as np
from astropy import constants, units
from astropy.cosmology import Planck18
from scipy import integrate

from gensed_base import gensed_base

M_sun = constants.M_sun.cgs.value
c = constants.c.cgs.value
pc = (1 * units.pc).to_value(units.cm)
G = constants.G.cgs.value
sigma_sb = constants.sigma_sb.cgs.value
h = constants.h.cgs.value
k_B = constants.k_B.cgs.value


class AGN:
    def __init__(self, t0: float, Mi: float, M_BH: float, lam: np.ndarray, edd_ratio: float, rng): # constructor
        self.lam = np.asarray(lam)
        self.t0 = t0
        self.rng = np.random.default_rng(rng)

        self.ME_dot = self.find_ME_dot(M_BH)
        self.MBH_dot = self.find_MBH_dot(self.ME_dot, M_BH, edd_ratio)
        self.Fnu_average = self.find_Fnu_average_standard_disk(self.MBH_dot, self.lam, M_BH)

        self.tau = self.find_tau_v(self.lam, Mi, M_BH)
        self.sf_inf = self.find_sf_inf(self.lam, Mi, M_BH)

        self.t = t0
        self.delta_m = self._random() * self.sf_inf

    def step(self, t):
        dt = t - self.t

        self.t = t

        self.delta_m = (
                self.delta_m * np.exp(-dt / self.tau)
                + self.sf_inf * np.sqrt(1 - np.exp(-2 * dt / self.tau)) * self._random()
        )

    def __call__(self, t):
        self.step(t)
        return self.Fnu

    @property
    def Fnu(self):
        return 10 ** (-0.4 * self.delta_m) * self.Fnu_average

    def _random(self):
        return self.rng.normal(size=self.lam.size)

    @staticmethod
    def find_Fnu_average(lam, M_BH, eddington_ratio):
        # Input wavelength as array.
        # Return baseline (average value).

        z = 0.2
        mu = Planck18.distmod(z).value
        F_av = 10 ** (-0.4 * (20 + 48.6 - mu))
        Fnu_ave = np.full_like(lam, F_av)
        return Fnu_ave
        # return 1e-29 * (lam / 5000e-8)**(-1/3)

    @staticmethod
    def find_ME_dot(M_BH):
    #return in g/s
        return 1.4e18*M_BH/M_sun

    @staticmethod
    def find_MBH_dot(ME_dot, M_BH, eddington_ratio):
        return ME_dot * eddington_ratio
    # don't need yita?

    @staticmethod
    def T_0(M, Mdot, r_in):
        return (2 ** (3 / 4) * (3 / 7) ** (7 / 4) * (G * M * Mdot / (np.pi * sigma_sb * r_in ** 3)) ** (1 / 4))

    @staticmethod
    def r_0(r_in):
        return ((7 / 6) ** 2 * r_in)

    @staticmethod
    def x_fun(nu, T0, r, r0):
        return (h * nu / (k_B * T0) * (r / r0) ** (3 / 4))

    def find_flux_standard_disk(self, Mdot, nu, rin, rout, i, d, M):
        T0 = self.T_0(M, Mdot, rin)
        r0 = self.r_0(rin)
        xin = self.x_fun(nu, T0, rin, r0)
        xout = self.x_fun(nu, T0, rout, r0)
        fun_integr = lambda x: (x ** (5 / 3)) / (np.exp(x) - 1)
    #     integ, inte_err = integrate.quad(fun_integr, xin, xout)
        integ, inte_err = integrate.quad(fun_integr, 0, np.inf)

        return ((16 * np.pi) / (3 * d ** 2) * np.cos(i) * (k_B * T0 / h) ** (8 / 3) * h * (nu ** (1 / 3)) / (c ** 2) * (
                    r0 ** 2) * integ)

    def find_Fnu_average_standard_disk(self, MBH_dot, lam, M_BH):

        flux_av = self.find_flux_standard_disk(MBH_dot, c/lam, rin=1, rout=1, i=0, d=10 * pc, M=M_BH) # do we need i, d?
        return flux_av


    @staticmethod
    def find_tau_v(lam, Mi=-23, M_BH=1e9 * M_sun):
        """Input frequency v in Hz, i band magnitude (default is -23), Black whole mass in g (defalt is 10^9 solar mass).
        Return timescale in s."""

        A = 2.4  # self.rng.normal(2.4, ...)
        B = 0.17
        C = 0.03
        D = 0.21
        # add C, D, BH_mass, Mi
        return 10 ** (A + B * np.log10(lam / (4000e-8))
                      + C * (Mi + 23) + D * np.log10(M_BH / (1e9 * M_sun)))  # e-8 angstrom

    @staticmethod
    def find_sf_inf(lam, Mi=-23, M_BH=1e9 * M_sun):
        """Input frequency in Hz, i band magnitude (default is -23), Black whole mass in g (defalt is 10^9 solar mass).

        Return Structure Function at infinity in mag."""

        A = -0.51
        B = -0.479
        C = 0.13
        D = 0.18

        return 10 ** (A + B * np.log10(lam / (4000e-8))
                      + C * (Mi + 23) + D * np.log10(M_BH / (1e9 * M_sun)))


class gensed_AGN(gensed_base):
    # round input rest time by 10**_trest_digits days
    _trest_digits = 8

    def __init__(self, PATH_VERSION, OPTMASK, ARGLIST, HOST_PARAM_NAMES):
        self.agn = None
        self.trest = None
        self.sed = None
        self.rng = np.random.default_rng(0)
        self.wavelen = 100
        self.wave = np.logspace(np.log10(100e-8), np.log10(20000e-8), self.wavelen)

    def _get_Flambda(self):
        return self.agn.Fnu * c / self.wave ** 2 * 1e-8

    def prepEvent(self, trest, external_id, hostparams):
        # trest is sorted
        self.trest = np.round(trest, self._trest_digits)
        self.agn = AGN(t0=self.trest[0], Mi=-23, M_BH=1e9 * M_sun, lam=self.wave, edd_ratio=0.1, rng=self.rng)
        self.sed = {self.trest[0]: self._get_Flambda()}
        # TODO: consider a case of repeated t, we usually have several t = 0
        for t in self.trest[1:]:
            self.agn.step(t)
            self.sed[t] = self._get_Flambda()

    def fetchSED_LAM(self):
        """
        Returns the wavelength vector
        """
        wave_aa = self.wave * 1e8
        # print('wave:',wave_aa)
        return wave_aa

    def fetchSED(self, trest, maxlam, external_id, new_event, hostparams):
        trest = round(trest, self._trest_digits)
        return self.sed[trest]

    def fetchParNames(self):
        return []

    def fetchParVals(self, varname):
        return 'SNANA'
