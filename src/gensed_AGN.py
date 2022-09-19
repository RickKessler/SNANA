import numpy as np
from astropy import constants
from astropy.cosmology import Planck18

M_sun = constants.M_sun.cgs.value
c = constants.c.cgs.value


class AGN:
    def __init__(self, t0: float, Mi: float, M_BH: float, lam: np.ndarray, rng):  # constructor
        self.lam = np.asarray(lam)
        self.t0 = t0
        self.rng = np.random.default_rng(rng)

        self.Fnu_average = self.find_Fnu_average(self.lam, M_BH, None)

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


class gensed_BAYESN:
    def __init__(self, PATH_VERSION, OPTMASK, ARGLIST, HOST_PARAM_NAMES):
        self.agn = None
        self.rng = np.random.default_rng(0)
        self.wavelen = 100
        self.wave = np.logspace(np.log10(100e-8), np.log10(20000e-8), self.wavelen)

    def fetchSED_NLAM(self):
        """
        Returns the length of the wavelength vector
        """
        return self.wavelen

    def fetchSED_LAM(self):
        """
        Returns the wavelength vector
        """
        wave_aa = self.wave * 1e8
        # print('wave:',wave_aa)
        return wave_aa.tolist()

    def fetchSED_BAYESN(self, trest, maxlam=5000, external_id=1, new_event=1, hostparams=''):
        if new_event:
            self.agn = AGN(t0=trest, Mi=-23, M_BH=1e9 * M_sun, lam=self.wave, rng=self.rng)
        else:
            self.agn.step(trest)
        Flambda = self.agn.Fnu * c / self.wave ** 2 * 1e-8
        return Flambda.tolist()

    def fetchParNames_BAYESN(self):
        return []

    def fetchNParNames_BAYESN(self):
        return 0

    def fetchParVals_BAYESN_4SNANA(self, varname):
        return 'SNANA'

