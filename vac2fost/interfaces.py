"""
Where interface classes live.

interfaces are the bulk of this package, and hold most capabilities,
they are used to translate (MPI-AMRVAC) output files to .fits input suited for MCFOST.

VtuFileInterface is the historic class, and will soon be deprecated
DatFileInterface is the future of this package (wip)
"""

# stdlib
import os
from pathlib import Path
from abc import ABC, abstractmethod

# non standard externals
import numpy as np
from astropy import units
from astropy.io import fits
from scipy.interpolate import interp1d, interp2d
import f90nml

# private externals
from vtk_vacreader import VacDataSorter

from .info import __version__
from .utils import shell_path
from .utils import IOinfo, DataInfo, GridShape
from .mcfost_utils import MINGRAINSIZE_µ, KNOWN_MCFOST_ARGS
from .mcfost_utils import get_mcfost_grid, write_mcfost_conf
from .logger import v2flogger as log



DEFAULT_UNITS = dict(distance2au=1.0, time2yr=1.0, mass2solar=1.0)

def read_amrvac_parfiles(parfiles: list, location: str = "") -> f90nml.Namelist:
    """Parse one, or a list of MPI-AMRVAC parfiles into a consistent
    configuration.

    <location> : pathlike of the directory where parfiles are found.
    Can be either a PathLike object or a str. The later can include
    "$" shell env variables such as "$HOME".

    This function replicates that of MPI-AMRVAC, with a patching logic:
    for parameters redundant across parfiles, only last values are kept,
    except "&filelist:base_filename", for which values are cumulated.
    """
    pathloc = shell_path(location)

    if isinstance(parfiles, (str, os.PathLike)):
        pfs = [parfiles]
    else:
        pfs = parfiles
    assert all([isinstance(pf, (str, os.PathLike)) for pf in pfs])

    confs = [f90nml.read((pathloc / pf).resolve()) for pf in pfs]
    conf_tot = f90nml.Namelist()
    for c in confs:
        conf_tot.patch(c)

    base_filename = "".join([c.get("filelist", {}).get("base_filename", "") for c in confs])
    assert base_filename != ""
    conf_tot["filelist"]["base_filename"] = base_filename
    return conf_tot



class AbstractInterface(ABC):
    """A data transforming class. Holds most functionalities useful to
    vac2fost.main()"""

    def __init__(self,
                 config_file: Path,
                 nums: int = None, # or any int-returning iterable
                 output_dir: Path = Path.cwd(),
                 dust_bin_mode: str = None,
                 read_gas_density=False,
                 read_gas_velocity=False,
                 settling=False,
                 axisymmetry=False):

        # input checking
        if not isinstance(config_file, (str, Path)):
            raise TypeError(config_file)
        if not isinstance(output_dir, (str, Path)):
            raise TypeError(output_dir)
        if axisymmetry and read_gas_velocity:
            raise NotImplementedError

        # attribute storage
        self._base_args = {
            'config_file': Path(config_file),
            'output_dir': Path(output_dir),
            'nums': nums,
            'dust_bin_mode': dust_bin_mode,
            'read_gas_density': read_gas_density,
            'read_gas_velocity': read_gas_velocity
        }

        self.read_gas_velocity = read_gas_velocity
        self.use_settling = settling
        self.use_axisymmetry = axisymmetry

        # parse configuration file
        self.conf = f90nml.read(config_file)
        if self.use_axisymmetry and self.conf['mcfost_output'].get("n_az", 2) > 1:
            log.warning("specified 'n_az'>1 but axisymmetry flag present, overriding n_az=1")
            self.conf["mcfost_output"].update({"n_az": 1})

        # init iteration counter
        if nums is None:
            nums = self.conf["amrvac_input"]["nums"]
        if isinstance(nums, int):
            nums = [nums]  # make it iterable
        nums = list(set(nums))  # filter out duplicates and sort them
        def _iter_nums():
            for n in nums:
                yield n
        self._iter_nums = _iter_nums()
        self.iter_count = 0
        self.iter_max = len(nums)
        self.current_num = next(self._iter_nums)

        hydro_data_dir = shell_path(self.conf["amrvac_input"]["hydro_data_dir"])
        if not hydro_data_dir.is_absolute():
            options = self.conf['amrvac_input']
            p1 = Path.cwd()
            p2 = (Path(config_file).parent/hydro_data_dir).resolve()

            if isinstance(options['config'], (list, tuple)):
                fi = options['config'][0]
            else:
                fi = options['config']

            found = [(p/fi).is_file() for p in (p1, p2)]
            if all(found) and p1 != p2:
                errmess = f'can not guess if path "{hydro_data_dir}" '
                errmess += "is relative to cwd or configuration file"
                raise FileNotFoundError(errmess)

            if not any(found):
                raise FileNotFoundError(hydro_data_dir/options['config'][0])

            p = (p1, p2)[found.index(True)]
            self.conf['amrvac_input'].update({'hydro_data_dir': p.resolve()})
        self.amrvac_conf = read_amrvac_parfiles(
            parfiles=self.conf['amrvac_input']['config'],
            location=self.conf['amrvac_input']['hydro_data_dir']
        )

        self._µsizes = None
        self._input_data = None
        self.output_grid = None
        self._new_2D_arrays = None
        self._new_3D_arrays = None
        self._rz_slice = None

        self._set_dbm(dust_bin_mode)

        if not self.io.OUT.directory.exists():
            os.makedirs(self.io.OUT.directory)
            log.warning(f"dir {self.io.OUT.directory} was created")

        if not self.conf.get("units"):
            log.warning(f"&units parameter list not found. Assuming {DEFAULT_UNITS}")
            self.conf["units"] = f90nml.Namelist(DEFAULT_UNITS)
        else:
            for k, v in DEFAULT_UNITS.items():
                if not self.conf["units"].get(k):
                    log.warning(f"&units:{k} parameter not found. Assuming default {v}")
                    self.conf["units"][k] = v

    # abstract bits
    @abstractmethod
    def load_input_data(self) -> None:
        """Set self._input_data"""

    @property
    @abstractmethod
    def g2d_ratio(self):
        """Parse the gas to dust ratio."""
        # todo: make this part of load_input_data()

    @property
    @abstractmethod
    def density_keys(self) -> list:
        """Get the ordered list of density keys (gas, ds1, ds2, ds3...)
        where 'ds' reads 'dust species'
        """
        # todo: make this part of load_input_data()

    @property
    @abstractmethod
    def io(self) -> IOinfo:
        pass


    # public facing methods
    def preroll_mcfost(self) -> None:
        mcfost_conf_file = self.io.OUT.directory / "mcfost_conf.para"

        if not mcfost_conf_file.is_file():
            # Create a complete mcfost conf file using (by decreasing priority)
            # - amrvac initial configuration : self._translate_amrvac_config()
            # - user specifications : self.conf['mcfost_output']
            # - defaults (defined in mcfost_utils.py)
            mcfost_parameters = {}
            mcfost_parameters.update(self._translate_amrvac_config())

            #Get unrecognized arguments found in mcfost_output
            unknown_args = []
            for arg in self.conf["mcfost_output"].keys():
                if not arg.lower() in KNOWN_MCFOST_ARGS:
                    unknown_args.append(arg)
            if unknown_args:
                raise ValueError(f'Unrecognized MCFOST argument(s): {unknown_args}')

            mcfost_parameters.update(self.conf['mcfost_output'])
            write_mcfost_conf(output_file=mcfost_conf_file,
                              custom_parameters=mcfost_parameters)
            log.info(f"successfully wrote {mcfost_conf_file}")

        target_grid = get_mcfost_grid(self, mcfost_conf_file,
                                      output_dir=self.io.OUT.directory,
                                      require_run=(self.iter_count == 0))
        self.output_grid = {
            "array": target_grid,
            # (nr, nphi) 2D grids
            "z-slice_r": target_grid[0, :, 0, :],
            # (nr, nz) 2D grid (z points do not depend on phi)
            "phi-slice_z": target_grid[1, 0, :, :],
            # vectors (1D arrays)
            "ticks_r": target_grid[0, 0, 0, :],
        }
        if target_grid.shape[0] > 2: # usually the case unless 2D axisym grid !
            self.output_grid.update({"z-slice_phi": target_grid[2, :, 0, :],
                                     "ticks_phi": target_grid[2, :, 0, 0]})

    def write_output(self) -> None:
        """Write a .fits file suited for MCFOST input."""
        dust_bin_selector = {
            "gas-only": np.zeros(1, dtype="int64"),
            "dust-only": 1 + self.grain_micron_sizes[1:].argsort(),
            "mixed": self.grain_micron_sizes.argsort()
        }[self._dust_binning_mode]

        if self.use_axisymmetry:
            gas_field = self._rz_slice[0]
            dust_fields = self._rz_slice[dust_bin_selector]
        else:
            gas_field = self.new_3D_arrays[0]
            dust_fields = self.new_3D_arrays[dust_bin_selector]

        suppl_hdus = []
        assert (len(dust_bin_selector) > 1) == (self._bin_dust)
        if len(dust_bin_selector) > 1:
            # mcfost requires an HDU with grain sizes only if more than one population is present
            suppl_hdus.append(
                fits.ImageHDU(self.grain_micron_sizes[dust_bin_selector])
            )

        header = {'read_n_a': 0} # automatic normalization of size-bins from mcfost param file.
        if self.read_gas_density:
            #devnote: add try statement here ?
            header.update(dict(gas_to_dust=self.amrvac_conf["usr_dust_list"]["gas2dust_ratio"]))
            suppl_hdus.append(fits.ImageHDU(gas_field))
            header.update(dict(read_gas_density=1))

        if self.read_gas_velocity:
            header.update(dict(read_gas_velocity=1))
            suppl_hdus.append(fits.ImageHDU(self.new_3D_gas_velocity))

        dust_densities_HDU = fits.PrimaryHDU(dust_fields)
        for k, v in header.items():
            # this is the canonical way to avoid HIERARCH-related warnings from astropy
            if len(k) > 8:
                k = f"HIERARCH {k}"
            dust_densities_HDU.header.append((k, v))

        with open(self.io.OUT.filepath, mode="wb") as fo:
            hdul = fits.HDUList(hdus=[dust_densities_HDU] + suppl_hdus)
            hdul.writeto(fo)
        log.info(f"successfully wrote {self.io.OUT.filepath}")

    def advance_iteration(self) -> None:
        """Reset output attributes and step to next output number."""
        self._new_2D_arrays = None
        self._new_3D_arrays = None
        self._rz_slice = None
        self.current_num = next(self._iter_nums)
        self.iter_count += 1

    # abusive properties ...
    @property
    def read_gas_density(self) -> bool:
        """Named after mcfost's option. Gas density is passed to mcfost only
        if required by user AND non-redundant.

        Clarification: if no gas density is passed, mcfost assumes
        that gas is traced by smallest grains. As "gas-only" and
        "mixed" modes make the same assumption, they would produce
        identical result without explicitly passing the gas density.
        """
        if not self._base_args["read_gas_density"]:
            rgd = False
        elif self._bin_gas:
            log.warning("specified read_gas_density but redundant"
                        f"with '{self._dust_binning_mode}' dbm, ignored")
            rgd = False
        else:
            rgd = True
        return rgd

    @property
    def grain_micron_sizes(self) -> np.ndarray:
        """Read grain sizes (assumed in [cm]), from AMRVAC parameters and
        convert to microns."""
        if self._µsizes is None:
            µm_sizes = np.empty(0)
            if self._bin_dust:
                µm_sizes = 1e4 * np.array(
                    self.amrvac_conf["usr_dust_list"]["grain_size_cm"])
                assert min(µm_sizes) > 0.1 #in case this triggers, review this code
            # always associate a grain size to the gas bin
            µm_sizes = np.insert(µm_sizes, 0, MINGRAINSIZE_µ)
            self._µsizes = µm_sizes
        return self._µsizes

    @property
    def aspect_ratio(self):
        """Dimensionless gas scale height implied by mcfost parameters"""
        mcfl = self.conf["mcfost_output"]
        return mcfl["scale_height"] / mcfl["reference_radius"]

    @property
    def input_grid(self) -> dict:
        """Store physical coordinates (vectors) about the input grid specifications."""
        ig = {
            "ticks_r": self._input_data.get_ticks("r") * self.conf["units"]["distance2au"],
            "ticks_phi": self._input_data.get_ticks("phi")
        }
        return ig

    # private methods
    def _set_dbm(self, dbm=None) -> None:
        """Define binning strategy
        - (gas-only)  : use only gas as a proxy for dust
        - (dust-only) : use only dust information
        - (mixed)     : use both, assuming gas traces the smallest grains
        """
        if dbm in ("dust-only", "gas-only", "mixed"):
            self._dust_binning_mode = dbm
        elif isinstance(dbm, str):
            raise ValueError(f"Unrecognized dbm value {dbm}")
        else: # automatic setting
            try: # todo: make this try block more robust
                smallest_gs_µm = 1e4* min(np.array(self.amrvac_conf['usr_dust_list']['grain_size_cm']))
            except KeyError:
                self._dust_binning_mode = "gas-only"
                reason = "could not find grain sizes"
            else:
                if smallest_gs_µm > MINGRAINSIZE_µ:
                    self._dust_binning_mode = "mixed"
                    reason = f"smallest size found > {MINGRAINSIZE_µ}µm"
            log.warning(f"switched to '{self._dust_binning_mode}' dbm ({reason})")

        # should dust fluids be passed to mcfost ?
        self._bin_dust = self._dust_binning_mode in {"dust-only", "mixed"}

        # should gas be passed to mcfost ?
        self._bin_gas = self._dust_binning_mode in {'gas-only', 'mixed'}

    def _estimate_dust_mass(self) -> float:
        """Estimate the total dust mass in the grid, in solar masses"""
        # devnote : this assumes a linearly spaced grid
        dphi = 2*np.pi / self.io.IN.gridshape.nphi
        rvect = self._input_data.get_ticks("r")
        dr = rvect[1] - rvect[0]
        cell_surfaces = dphi/2 * ((rvect + dr/2)**2 - (rvect - dr/2)**2)

        if self._dust_binning_mode == "gas-only":
            keys = ["rho"]
        else:
            keys = [k for k, _ in self._input_data if "rhod" in k]
        mass = 0.0
        for key in keys:
            mass += np.sum([cell_surfaces * self._input_data[key][:, i]
                            for i in range(self.io.IN.gridshape.nphi)])
        if self._dust_binning_mode == "gas-only":
            mass /= self.g2d_ratio
        mass *= self.conf["units"]["mass2solar"]
        return mass

    def _translate_amrvac_config(self) -> dict:
        parameters = {}

        # Zone
        mesh = self.amrvac_conf['meshlist']
        conv2au = self.conf["units"]["distance2au"]
        parameters.update({
            'rin': mesh['xprobmin1']*conv2au,
            'rout': mesh['xprobmax1']*conv2au,
            'maps_size': 2*mesh['xprobmax1']*conv2au,
        })

        if self._bin_dust:
            parameters.update({
                "gas_to_dust_ratio": self.g2d_ratio,
                "dust_mass": self._estimate_dust_mass()
            })
            # Grains
            sizes_µm = self.grain_micron_sizes
            parameters.update({
                # min/max grain sizes in microns
                'sp_min': min(1e-1, min(sizes_µm)),
                'sp_max': max(1e3, max(sizes_µm)),
            })
        #Star
        try:
            parameters.update({"star_mass": self.amrvac_conf["disk_list"]["central_mass"]})
        except KeyError:
            log.warning("&disk_list not found. Assuming default values")
        return parameters



    def _interpolate2D(self, datakey: str) -> np.ndarray:
        """Transform a polar field from MPI-AMRVAC coords to mcfost coords"""
        interpolator = interp2d(
            self.input_grid["ticks_phi"],
            self.input_grid["ticks_r"],
            self._input_data[datakey],
            kind="cubic",
            copy=False # test
        )
        return interpolator(self.output_grid["ticks_phi"], self.output_grid["ticks_r"])

    def _interpolate1D(self, datakey: str) -> np.ndarray:
        interpolator = interp1d(
            self.input_grid["ticks_r"],
            self._input_data[datakey][:, 0], # radial profile
            kind="cubic",
            copy=False,
            fill_value="extrapolate"
        )
        return interpolator(self.output_grid["ticks_r"])


    def gen_rz_slice(self) -> None:
        radial_profiles = np.array([self._interpolate1D(datakey=k) for k in self.density_keys])
        oshape = self.io.OUT.gridshape
        nr, nz = oshape.nr, oshape.nz

        nbins = len(radial_profiles)
        self._rz_slice = np.zeros((nbins, nz, nr))
        for ir, r in enumerate(self.output_grid["ticks_r"]):
            z_vect = self.output_grid["phi-slice_z"][:, ir]
            gas_height = r * self.aspect_ratio
            for i_bin, (grain_µsize, rprofile) in enumerate(zip(self.grain_micron_sizes,
                                                                radial_profiles)):
                H = gas_height
                if self.use_settling:
                    H *= (grain_µsize / MINGRAINSIZE_µ)**(-0.5)
                gaussian = np.exp(-z_vect**2/ (2*H**2)) / (np.sqrt(2*np.pi) * H)
                self._rz_slice[i_bin, :, ir] = gaussian * rprofile[ir]

    def gen_2D_arrays(self) -> None:
        """Interpolate input data density fields from input coords to output coords"""
        self._new_2D_arrays = np.array([self._interpolate2D(datakey=k) for k in self.density_keys])
        assert self._new_2D_arrays[0].shape == (self.io.OUT.gridshape.nr,
                                                self.io.OUT.gridshape.nphi)

    def gen_3D_arrays(self) -> None:
        """Interpolate input data onto full 3D output grid"""
        oshape = self.io.OUT.gridshape
        nr, nphi, nz = oshape.nr, oshape.nphi, oshape.nz

        nbins = len(self.new_2D_arrays)
        self._new_3D_arrays = np.zeros((nbins, nphi, nz, nr))
        for ir, r in enumerate(self.output_grid["ticks_r"]):
            z_vect = self.output_grid["phi-slice_z"][nz:, ir].reshape(1, nz)
            gas_height = r * self.aspect_ratio
            for i_bin, grain_µsize in enumerate(self.grain_micron_sizes):
                surface_density = self.new_2D_arrays[i_bin, ir, :]
                H = gas_height
                if self.use_settling:
                    H *= (grain_µsize / MINGRAINSIZE_µ)**(-0.5)
                gaussian = np.exp(-z_vect**2/ (2*H**2)) / (np.sqrt(2*np.pi) * H)
                #todo: numpy ellipsis ? "..."
                self._new_3D_arrays[i_bin, :, :, ir] = \
                    gaussian * surface_density.reshape(nphi, 1)

    @property
    def new_2D_arrays(self) -> list:
        #todo: rename me
        '''Last minute generation is used if required'''
        if self._new_2D_arrays is None:
            self.gen_2D_arrays()
        return self._new_2D_arrays

    @property
    def new_3D_arrays(self) -> list:
        '''Last minute generation is used if required'''
        #todo: rename me
        if self._new_3D_arrays is None:
            self.gen_3D_arrays()
        return self._new_3D_arrays

    @property
    def new_3D_gas_velocity(self) -> np.ndarray:
        """Derive the 3D velocity field for gas velocity, in km/s"""
        rho, mr, mphi = map(self._interpolate2D, ["rho", "m1", "m2"])
        vr, vphi = map(lambda x: x/rho, [mr, mphi])
        phig = self.output_grid["z-slice_phi"].transpose()
        vx = vr * np.cos(phig) - vphi * np.sin(phig)
        vy = vr * np.sin(phig) + vphi * np.cos(phig)

        # transform to 3D
        nz = self.io.OUT.gridshape.nz
        vx, vy = map(lambda a: np.stack([a]*nz, axis=1), [vx, vy])
        vz = np.zeros(vx.shape)
        oshape = self.io.OUT.gridshape
        for v in (vx, vy, vz):
            np.testing.assert_array_equal(v.shape, (oshape.nr, oshape.nz, oshape.nphi))

        # unit conversion
        conv = self.conf["units"]
        dimvel = conv["distance2au"]*units.au / (conv["time2yr"]*units.yr)
        vel2kms = dimvel.to(units.m / units.s).value
        velarr = np.stack([vx, vy, vz], axis=3) * vel2kms
        return velarr.transpose()


class VtuFileInterface(AbstractInterface):
    """An interface dedicated to fixed-resolution vtu files (to be deprecated !)"""
    @property
    def io(self) -> IOinfo:
        """Give up-to-date information on data location and naming (.i: input, .o: output)"""
        vtufile_name = "".join([self.amrvac_conf["filelist"]["base_filename"],
                                str(self.current_num).zfill(4),
                                ".vtu"])

        geomdefs = {"nr": 1, "nphi": 2}
        _input = DataInfo(
            directory=shell_path(self.conf["amrvac_input"]["hydro_data_dir"]).resolve(),
            filename=vtufile_name,
            gridshape=GridShape(**{k: self.amrvac_conf["meshlist"][f"domain_nx{n}"]
                                   for k, n in geomdefs.items()})
        )

        trad_keys = {"nr": "n_rad", "nphi": "n_az", "nz": "nz"}
        _output = DataInfo(
            directory=Path(self._base_args["output_dir"]),
            filename=_input.filestem+".fits",
            gridshape=GridShape(**{k1: self.conf["mcfost_output"][k2]
                                   for k1, k2 in trad_keys.items()})
        )
        return IOinfo(IN=_input, OUT=_output)

    def load_input_data(self) -> None:
        """Use vtkvacreader.VacDataSorter to load AMRVAC data"""
        self._input_data = VacDataSorter(
            file_name=str(self.io.IN.filepath),
            shape=(self.io.IN.gridshape.nr, self.io.IN.gridshape.nphi)
        )
        log.info(f"successfully loaded {self.io.IN.filepath}")

    @property
    def density_keys(self) -> list:
        return sorted(filter(lambda k: "rho" in k, self._input_data.fields.keys()))

    @property
    def g2d_ratio(self):
        """Gas to dust ratio"""
        return self.amrvac_conf["usr_dust_list"].get("gas2dust_ratio", 100.)


class DatFileInterface(AbstractInterface):
    @property
    def io(self) -> IOinfo:
        """Give up-to-date information on data location and naming (.i: input, .o: output)"""
        pass

    def load_input_data(self) -> None:
        """wip"""
        import yt
        from yt import mylog as ytlogger
        ytlogger.setLevel(50) # silence all but critical messages from yt
        dims = self.io.IN.gridshape
        print(dims)
        ag = ds.arbitrary_grid(ds.domain_left_edge, ds.domain_right_edge, dims=dims)
        d = ag["rho"].squeeze().to_ndarray()
        # TODO: use yt to set self._input_data here
        #log.info(f"successfully loaded {filename}")
