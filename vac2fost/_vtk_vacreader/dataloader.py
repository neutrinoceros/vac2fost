from pathlib import Path
from collections import defaultdict

import numpy as np
import f90nml

from vtk_vacreader import VacDataSorter, AugmentedVacDataSorter


class VacDataLoader:
    """A container class that allows reuse of <basename> and <shape> arguments.
    Load data arrays as usage needed.
    """

    _basedatasorter = VacDataSorter

    def __init__(self, basename: str, shape: tuple = None, sim_params: f90nml.Namelist = None):
        self.basename = basename
        self.shape = shape
        self.sim_params = sim_params
        self._loaded = defaultdict(bool)
        self._data = {}

    def __getitem__(self, key: int):
        targetfile = f"{self.basename}{str(key).zfill(4)}.vtu"
        if not self._loaded[key]:
            if not Path(targetfile).exists():
                raise FileNotFoundError(targetfile)
            else:
                self._data.update(
                    {
                        key: self.__class__._basedatasorter(
                            sim_params=self.sim_params, file_name=targetfile, shape=self.shape
                        )
                    }
                )
                self._loaded[key] = True
        return self._data[key]

    @property
    def meshgrid(self):
        # devnote : this will fail if the initial snapshot is not present
        # also, it implies that we load a snapshot just to extract the meshgrid... not very efficient
        if self._meshgrid is None:
            self._meshgrid = self[0].get_meshgrid()
        return self._meshgrid

    @property
    def phig(self) -> np.ndarray:
        return self.meshgrid[0]

    @property
    def phig_deg(self) -> np.ndarray:
        return np.rad2deg(self.meshgrid[0])

    @property
    def phig_deg_center(self) -> np.ndarray:
        return self.phig_deg - 180

    @property
    def rg(self) -> np.ndarray:
        return self.meshgrid[1]

    def get_phytime(self, num: int) -> float:
        """convert output number into physical time (code units)"""
        return num * self.sim_params["savelist"]["dtsave_dat"]


class AugmentedVacDataLoader(VacDataLoader):
    """Add a few functionnalities on top of VacDataLoader that are not
    general enough to include in the parent class"""

    _basedatasorter = AugmentedVacDataSorter

    def __init__(self, **args):
        super(__class__, self).__init__(**args)
        self._meshgrid = None
        self._vK = None

    @property
    def keplerian_velocity(self) -> np.ndarray:
        if self._vK is None:
            rv = self[0].get_ticks(axis=0)
            M_s = self.sim_params["disk_list"]["central_mass"]
            G = 4 * np.pi ** 2 * self.sim_params["disk_list"]["ref_radius"] ** 3 / M_s
            self._vK = np.sqrt(G * M_s) * rv ** -0.5
        return self._vK

    @property
    def keplerian_velocity_grid(self) -> np.ndarray:
        return np.stack([self.keplerian_velocity for i in range(self.shape[1])], axis=1)
