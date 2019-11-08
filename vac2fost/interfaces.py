"""
Where interface classes live.

interfaces are the bulk of this package, and hold most capabilities,
they are used to translate (MPI-AMRVAC) output files to .fits input suited for MCFOST.
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
from .mcfost_utils import MINGRAINSIZE_mum, KNOWN_MCFOST_ARGS
from .mcfost_utils import get_mcfost_grid, write_mcfost_conf
from .logger import v2flogger as log



DEFAULT_UNITS = dict(distance2au=1.0, time2yr=1.0, mass2solar=1.0)

def read_amrvac_parfiles(parfiles: list, location: str = "") -> f90nml.Namelist:
    """Parse one, or a list of MPI-AMRVAC parfiles into a consistent configuration.

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

    def __init__(self, conf_file: Path, override: dict = None, output_dir: Path = None):
        # python 3.8: make conf_file positional only

        # input checking
        if not isinstance(conf_file, (str, Path)):
            raise TypeError(conf_file)
        if output_dir is None:
            output_dir = Path.cwd()
            log.warning("no output_dir provided, outputs will be written to current work directory")
        elif not isinstance(output_dir, (str, Path)):
            raise TypeError(output_dir)
        if override is None:
            override = {}
        elif not isinstance(override, dict):
            raise TypeError(override)

        # parse configuration
        self._output_dir = Path(output_dir)
        self._output_conf_file = self._output_dir / "vac2fost.nml.backup"
        if self._output_conf_file == Path(conf_file):
            err = f"{self._output_conf_file} is a reserved file name for vac2fost to output. "
            err += "It can not be used as input."
            raise RuntimeError(err)

        self.conf_file = Path(conf_file)
        self.conf = f90nml.read(conf_file)
        self.conf.patch(override)

        flags = self.conf.get("flags", {})
        self._use_settling = flags.get("settling", False)
        self._use_axisymmetry = flags.get("axisymmetry", False)
        self._read_gas_velocity = flags.get("read_gas_velocity", False)

        # handle special cases
        if self._use_axisymmetry and self._read_gas_velocity:
            raise NotImplementedError
        if self._use_axisymmetry and self.conf["mcfost_output"].get("n_az", 2) > 1:
            log.warning("specified n_az > 1 but axisymmetry flag present, overriding n_az = 1")
            self.conf["mcfost_output"].update({"n_az": 1})

        # init iteration counter
        nums = self.conf["amrvac_input"]["nums"] # mandatory argument
        if isinstance(nums, int):
            nums = [nums]  # make it iterable
        nums = list(set(nums))  # filter out duplicates and sort them
        def _iter_nums():
            for n in nums:
                yield n
        self._iter_nums = _iter_nums()
        self._iter_count = 0
        self._iter_max = len(nums)
        self.current_num = next(self._iter_nums)

        hydro_data_dir = shell_path(self.conf["amrvac_input"]["hydro_data_dir"])
        if not hydro_data_dir.is_absolute():
            options = self.conf['amrvac_input']
            p1 = Path.cwd()
            p2 = (Path(conf_file).parent/hydro_data_dir).resolve()

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
            log.warning("Relative path found for hydro_data_dir, overriding to absolute path.")
            self.conf['amrvac_input'].update({'hydro_data_dir': str(p.resolve())})
        self.amrvac_conf = read_amrvac_parfiles(
            parfiles=self.conf['amrvac_input']['config'],
            location=self.conf['amrvac_input']['hydro_data_dir']
        )

        self._mumsizes = None
        self._input_data = None
        self._density_keys = None
        self.output_grid = None

        self._parse_dust_properties()

        read_gas_density = flags.get("read_gas_density", False)
        if read_gas_density and self._bin_gas:
            log.warning(f"Found redundancy: with dust_bin_mode='{self._dust_bin_mode}'"
                        "flag read_gas_density will be ignored.")
            self._read_gas_density = False
        else:
            # Clarification: if no gas density is passed, mcfost assumes
            # that gas is traced by smallest grains. As "gas-only" and
            # "mixed" modes make the same assumption, they would produce
            # identical result without explicitly passing the gas density.
            self._read_gas_density = read_gas_density

        if not self._output_dir.exists(): # python 3.8: :=
            os.makedirs(self._output_dir)
            log.warning(f"dir {self._output_dir} was created")

        # sanitize conf
        default_mcfost_output = {"scale_height": 0.05, "reference_radius": 1}
        for namelist, defaults in (("units", DEFAULT_UNITS),
                                   ("mcfost_output", default_mcfost_output)):
            if namelist not in self.conf:
                log.warning(f"&{namelist} parameter list not found. Overriding {defaults}")
                self.conf[namelist] = f90nml.Namelist(defaults)
            else:
                for k, v in defaults.items():
                    if not self.conf[namelist].get(k):
                        log.warning(f"&{namelist}:{k} parameter not found. Overriding default {v}")
                        self.conf[namelist][k] = v

    # abstract bits
    @abstractmethod
    def load_input_data(self) -> None:
        """Set self._input_data"""
        # this needs to set the following attributes
        #self._input_data
        #self._density_keys

    @property
    @abstractmethod
    def io(self) -> IOinfo:
        """Store input/ouput directories and grids information"""

    @property
    @abstractmethod
    def input_grid(self) -> dict:
        """Store physical coordinates (vectors) about the input grid specifications."""

    # public methods, for direct usage in vac2fost.main()
    def preroll_mcfost(self) -> None:
        """Output mcfost parafile and store the grid"""
        mcfost_conf_file = self.io.OUT.directory / "mcfost_conf.para"

        if not mcfost_conf_file.is_file():
            # Create a complete mcfost conf file using (by decreasing priority)
            # - amrvac initial configuration : self._translate_amrvac_config()
            # - user specifications : self.conf['mcfost_output']
            # - defaults (defined in mcfost_utils.py)
            mcfost_parameters = {}
            mcfost_parameters.update(self._translate_amrvac_config())
            mcfost_parameters.update(self.conf['mcfost_output'])

            #Star mass is a special case
            mstar = self.conf["mcfost_output"].get("mstar", None)
            if mstar is None:
                log.warning("&mcfost_output: Mstar not found. Assuming default value.")
                mstar = 1.0
            elif isinstance(mstar, str):
                namelist, param = mstar.split(".")
                mstar = self.amrvac_conf[namelist][param] * self.conf["units"]["mass2solar"]
            mcfost_parameters.update({"mstar": mstar})

            #Get unrecognized arguments found in mcfost_output
            unknown_args = []
            for arg in self.conf["mcfost_output"].keys():
                if not arg.lower() in KNOWN_MCFOST_ARGS:
                    unknown_args.append(arg)
            if unknown_args:
                raise ValueError(f'Unrecognized MCFOST argument(s): {unknown_args}')

            write_mcfost_conf(output_file=mcfost_conf_file,
                              custom_parameters=mcfost_parameters)
            log.info(f"successfully wrote {mcfost_conf_file}")

        grid = get_mcfost_grid(mcfost_conf_file,
                               output_dir=self.io.OUT.directory,
                               require_run=(self._iter_count == 0))

        self.output_grid = {
            "array": grid,
            # r coords
            "ticks_r": grid[0, 0, 0, :],
            "z-slice_r": grid[0, :, 0, :],
            "phi-slice_r": grid[0, 0, :, :],
            # z coords
            # note: we purposedly do not define a "ticks-z" 1D array because its value varies with r
            "phi-slice_z": grid[1, 0, :, :]
        }
        if grid.shape[0] > 2: # usually the case unless 2D axisym grid !
            self.output_grid.update({"ticks_phi": grid[2, :, 0, 0],
                                     "z-slice_phi": grid[2, :, 0, :]})

    def write_output(self) -> None:
        """Write a .fits file suited for MCFOST input."""
        dust_bin_selector = {
            "gas-only": np.zeros(1, dtype="int64"),
            "dust-only": 1 + self._grain_micron_sizes[1:].argsort(),
            "mixed": self._grain_micron_sizes.argsort()
        }[self._dust_bin_mode]

        output_ndarray = self.get_output_ndarray()
        gas_field = output_ndarray[0]
        dust_fields = output_ndarray[dust_bin_selector]

        suppl_hdus = []
        assert (len(dust_bin_selector) > 1) == (self._bin_dust)
        if len(dust_bin_selector) > 1:
            # mcfost requires an HDU with grain sizes only if more than one population is present
            suppl_hdus.append(
                fits.ImageHDU(self._grain_micron_sizes[dust_bin_selector])
            )

        header = {'read_n_a': 0} # automatic normalization of size-bins from mcfost param file.
        if self._read_gas_density:
            header.update(dict(gas_to_dust=self._gas_to_dust_ratio))
            suppl_hdus.append(fits.ImageHDU(gas_field))
            header.update(dict(read_gas_density=1))

        if self._read_gas_velocity:
            vel_array = self.get_gas_velocity_ndarray()
            header.update(dict(read_gas_velocity=1))
            suppl_hdus.append(fits.ImageHDU(vel_array))

        dust_densities_HDU = fits.PrimaryHDU(dust_fields)
        for k, v in header.items():
            # this is the canonical way to avoid HIERARCH-related warnings from astropy
            if len(k) > 8:
                k = f"HIERARCH {k}"
            dust_densities_HDU.header.append((k, v))

        if self._iter_count == 0:
            self.conf.write(self._output_conf_file, force=True)
            with open(self._output_conf_file, mode="at") as stream:
                stream.write(f"! automatically generated with vac2fost {__version__}\n")
                stream.write("! this file is self-contained and can be used for reproduction\n")
                stream.write("! WARNING: rename this file before running vac2fost with it")
            log.info(f"wrote {self._output_conf_file.resolve()}")

        with open(self.io.OUT.filepath, mode="wb") as fo:
            hdul = fits.HDUList(hdus=[dust_densities_HDU] + suppl_hdus)
            hdul.writeto(fo)
        log.info(f"wrote {self.io.OUT.filepath}")

    def advance_iteration(self) -> None:
        """Step to next output number."""
        self.current_num = next(self._iter_nums)
        self._iter_count += 1

    @property
    def iter_frac(self):
        """Visual hint for iteration advancement."""
        return f"{self._iter_count+1}/{self._iter_max}"

    @property
    def iter_last(self):
        """Whether or not we're at last iteration."""
        return self._iter_count == self._iter_max



    # private methods
    def _parse_dust_properties(self) -> None:
        """
        Parse grain sizes and define binning strategy
        - (gas-only)  : use only gas as a proxy for dust
        - (dust-only) : use only dust information
        - (mixed)     : use both, assuming gas traces the smallest grains
        """
        try:
            dbm = self.conf["flags"]["dust_bin_mode"]
        except KeyError:
            dbm = None
        if dbm is not None and dbm not in ("dust-only", "mixed", "gas-only"):
            raise ValueError(f"Unrecognized dbm value {dbm}")

        mum_sizes = np.empty(0)
        auto_dbm = None

        # automatic setup
        if "dust" in self.conf:
            grain_size2micron = self.conf["dust"].get("grain_size2micron", 1.0)
            grain_sizes = self.conf["dust"]["grain_sizes"]
            if isinstance(grain_sizes, str):
                namelist, item = grain_sizes.split(".")
                grain_sizes = self.amrvac_conf[namelist][item]
            grain_sizes = np.array(grain_sizes)
            mum_sizes = grain_size2micron * grain_sizes

        if dbm is not None:
            self._dust_bin_mode = dbm
        else:
            if len(mum_sizes) < 1:
                auto_dbm = "gas-only"
                reason = "dust parameters not found"
            elif min(mum_sizes) > MINGRAINSIZE_mum:
                auto_dbm = "mixed"
                reason = f"smallest grain size > threshold ({MINGRAINSIZE_mum} µm)"
            else:
                auto_dbm = "dust-only"
                reason = f"smallest grain size < threshold ({MINGRAINSIZE_mum} µm)"
            self._dust_bin_mode = auto_dbm
            log.warning(f"switched to '{self._dust_bin_mode}' dbm ({reason})")

        if self._dust_bin_mode == "gas-only":
            mum_sizes = np.empty(0)

        # always associate a grain size to the gas bin
        self._grain_micron_sizes = np.insert(mum_sizes, 0, MINGRAINSIZE_mum)

        # should dust fluids be passed to mcfost ?
        self._bin_dust = self._dust_bin_mode in {"dust-only", "mixed"}

        # should gas be passed to mcfost ?
        self._bin_gas = self._dust_bin_mode in {'gas-only', 'mixed'}

        if self._bin_dust:
            # parse gas to dust mass ratio
            dustlist = self.conf["dust"]
            if "gas_to_dust_ratio" and "dust_to_gas_ratio" in dustlist:
                raise RuntimeError("Can not set both 'gas_to_dust_ratio' and 'dust_to_gas_ratio'")
            if "gas_to_dust_ratio" in dustlist:
                g2d = dustlist["gas_to_dust_ratio"]
            elif "dust_to_gas_ratio" in dustlist:
                g2d = 1/dustlist["dust_to_gas_ratio"]
            else:
                g2d = 100
                log.warning("Could not find 'gas_to_dust_ratio', defaulting to 100")
            self._gas_to_dust_ratio = g2d
        else:
            self._gas_to_dust_ratio = None

    def _estimate_dust_mass(self) -> float:
        """Estimate the total dust mass in the grid, in solar masses"""
        # devnote : this assumes a linearly spaced grid
        dphi = 2*np.pi / self.io.IN.gridshape.nphi
        rvect = self._input_data.get_ticks("r")
        dr = rvect[1] - rvect[0]
        cell_surfaces = dphi/2 * ((rvect + dr/2)**2 - (rvect - dr/2)**2)

        if self._dust_bin_mode == "gas-only":
            keys = ["rho"]
        else:
            keys = [k for k, _ in self._input_data if "rhod" in k]
        mass = 0.0
        for key in keys:
            mass += np.sum([cell_surfaces * self._input_data[key][:, i]
                            for i in range(self.io.IN.gridshape.nphi)])
        if self._dust_bin_mode == "gas-only":
            mass /= 100
        mass *= self.conf["units"]["mass2solar"]
        return mass

    def _translate_amrvac_config(self) -> dict:
        """Get some mcfost parameters directly from amrvac."""
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
                "gas_to_dust_ratio": self._gas_to_dust_ratio,
                "dust_mass": self._estimate_dust_mass()
            })
            # Grains
            parameters.update({
                # min/max grain sizes in microns
                'sp_min': min(1e-1, min(self._grain_micron_sizes)),
                'sp_max': max(1e3, max(self._grain_micron_sizes)),
            })
        return parameters


    # output generation
    def get_output_ndarray(self) -> np.ndarray:
        """Compute the approriate density array with shape (nbins, nphi, nz, nr),
        to be later written as a HDUPrimary
        """
        nbins = len(self._density_keys)
        oshape = self.io.OUT.gridshape
        nr, nphi, nz = oshape.nr, oshape.nphi, oshape.nz

        r_profile_densities = np.zeros((nbins, nr))
        phi_slice_densities = np.zeros((nbins, 1, nz, nr))

        new_plane_densities = np.zeros((nbins, nr, nphi))
        full3D_densities = np.zeros((nbins, nphi, nz, nr))

        if self._use_axisymmetry:
            interpolate = self._interpolate1D
            # those are references, not copies
            hyperplane_densities = r_profile_densities
            output_ndarray = phi_slice_densities

        else:
            interpolate = self._interpolate2D
            # those are references, not copies
            hyperplane_densities = new_plane_densities
            output_ndarray = full3D_densities
        hyperplane_densities[:] = np.array([interpolate(datakey=k) for k in self._density_keys])


        #dimensionless gas scale height implied by mcfost parameters
        aspect_ratio = self.conf["mcfost_output"]["scale_height"] \
                       / self.conf["mcfost_output"]["reference_radius"]

        for ir, r in enumerate(self.output_grid["ticks_r"]):
            if self._use_axisymmetry: # devnote: unify these two lines ?
                z_vect = self.output_grid["phi-slice_z"][:, ir]
            else:
                z_vect = self.output_grid["phi-slice_z"][nz:, ir]
            gas_height = r * aspect_ratio
            for ibin, grain_mumsize in enumerate(self._grain_micron_sizes):
                hpd = hyperplane_densities[ibin, ir, ...]
                H = gas_height
                if self._use_settling:
                    H *= (grain_mumsize / MINGRAINSIZE_mum)**(-0.5)
                gaussian = np.exp(-z_vect**2/ (2*H**2)) / (np.sqrt(2*np.pi) * H)
                output_ndarray[ibin, ..., ir] = np.outer(hpd, gaussian)
        return output_ndarray

    def get_gas_velocity_ndarray(self) -> np.ndarray:
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

    def _interpolate2D(self, datakey: str) -> np.ndarray:
        """Transform a polar field (r-phi) from amrvac to mcfost coords."""
        interpolator = interp2d(
            self.input_grid["ticks_phi"],
            self.input_grid["ticks_r"],
            self._input_data[datakey],
            kind="cubic",
            copy=False
        )
        return interpolator(self.output_grid["ticks_phi"], self.output_grid["ticks_r"])

    def _interpolate1D(self, datakey: str) -> np.ndarray:
        """Transform a polar slice (r profile) from amrvac to mcfost coords."""
        interpolator = interp1d(
            self.input_grid["ticks_r"],
            self._input_data[datakey][:, 0], # radial profile
            kind="cubic",
            copy=False,
            fill_value="extrapolate"
        )
        return interpolator(self.output_grid["ticks_r"])



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
            directory=Path(self._output_dir),
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
        # ordered list of density keys (gas, ds1, ds2, ds3...)
        # where 'ds' reads 'dust species'
        self._density_keys = sorted(filter(lambda k: "rho" in k, self._input_data.fields.keys()))
        log.info(f"successfully loaded {self.io.IN.filepath}")

    @property
    def input_grid(self) -> dict:
        """Describe the amrvac grid."""
        ig = {"ticks_r": self._input_data.get_ticks("r") * self.conf["units"]["distance2au"],
              "ticks_phi": self._input_data.get_ticks("phi")}
        return ig