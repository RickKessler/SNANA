from abc import ABC, abstractmethod
from typing import Sequence, Tuple

import numpy as np
import numpy.typing as npt


class gensed_base(ABC):
    @abstractmethod
    def __init__(self, PATH_VERSION: str, OPTMASK: int, ARGLIST: str, HOST_PARAM_NAMES: str):
        pass


    @abstractmethod
    def fetchSED_LAM(self) -> npt.ArrayLike:
        """
        Numpy array of dtype float64 with wavelengths in Angstrom
        """
        raise NotImplementedError

    def _fetchSED_LAM(self) -> np.ndarray:
        """Wrapper of fetchSED_LAM to call from C"""
        return np.asarray(self.fetchSED_LAM(), dtype=np.float64)

    @abstractmethod
    def fetchSED(self, trest: float, maxlam: int, external_id: int, new_event: int, hostpars: Tuple[float]) -> npt.ArrayLike:
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
        hostpars : tuple[float]
             Host parameters corresponded to HOST_PARAM_NAMES

        Returns
        -------
        A numpy array of dtype float64 the same length as fetchSED_LAM containing the flux
        observed from 10 pc in erg / s / cm^2 / Angstrom
        """
        raise NotImplementedError

    def _fetchSED(self, *args, **kwargs) -> np.ndarray:
        """Wrapper of fetchSED to call from C"""
        return np.asarray(self.fetchSED(*args, **kwargs), dtype=np.float64)

    @abstractmethod
    def fetchParNames(self) -> Sequence[str]:
        """
        Returns the names of model parameters
        """
        raise NotImplementedError

    @abstractmethod
    def fetchParVals(self, varname: str) -> float:
        """
        Returns the value of parameter 'varname'

        Parameters
        ----------
        varname : str
             A parameter name from self.parameter_names
        """
        raise NotImplementedError

