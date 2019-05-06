#!/usr/bin/env python3
"""A conversion facility for MPI-AMRVAC (.vtu) to MCFOST (.fits)

This is a Python package (from vac2fost import main as vac2fost), and
also a command line script.  Run `python vac2fost.py --help` for
documentation on command line usage.

The main algorithm runs the following steps
  a) load AMRVAC data with vtk_vacreader.VacDataSorter(), sort it as 2D arrays
  b) dry-run MCFOST to get the exact output
  c) interpolate data to the target grid
  d) convert to 3D (gaussian redistribution of density)
  e) collect, sort and write output data to a fits file


Known limitations
  1) AMR grids are not supported (.vtu files need to be converted to uniform grids)
  2) .vtu are assumed to use polar coordinates (r-phi 2D)
  3) 2D interpolation does not account for the curvature of polar cells


Disclaimer
  This package is using Python3.7 syntax/features and will not be made backward
  compatible with older versions of Python.
"""



# Imports
# =======================================================================================
# stdlib
import os
import sys
from pathlib import Path
from argparse import ArgumentParser

# non standard externals
import numpy as np
from astropy import units
from astropy.io import fits
from scipy.interpolate import interp2d
import f90nml

# private externals
from vtk_vacreader import VacDataSorter
from vac2fost.info import __version__
from vac2fost.utils import colorama, RED, CYAN, BOLD
from vac2fost.utils import shell_path, wait_for_ok, get_prompt_size, decorated_centered_message
from vac2fost.utils import IOinfo, DataInfo, GridShape
from vac2fost.mcfost_utils import MINGRAINSIZE_µ, KNOWN_MCFOST_ARGS
from vac2fost.mcfost_utils import get_mcfost_grid, write_mcfost_conf




# Globals ===============================================================================
DEFAULT_UNITS = dict(distance2au=1.0, time2yr=1.0, mass2solar=1.0)



# Defintions ============================================================================
def generate_conf_template() -> f90nml.Namelist:
    """Generate a template namelist object with comments instead of default values"""
    amrvac_list = dict(
        hydro_data_dir="path/to/output/data/directory",
        config="relative/to/<hydro_data_dir>/path/to/amrvac/config/file[s]",
        nums=0
    )

    mcfost_list = dict(
        n_rad=128, n_rad_in=4, n_az=128, nz=10,
        # aspect ratio is implied by those parameters
        flaring_exp=1.125,
        reference_radius=100.0,  # [a.u.]
        scale_height=1.0,  # [a.u.], at defined at reference_radius
    )
    sublists = {
        "amrvac_input": amrvac_list,
        "units": DEFAULT_UNITS,
        "mcfost_output": mcfost_list
    }
    template = f90nml.Namelist({k: f90nml.Namelist(v) for k, v in sublists.items()})
    return template

def read_amrvac_parfiles(parfiles: list, location: str = "") -> f90nml.Namelist:
    """Parse one, or a list of MPI-AMRVAC parfiles into a consistent
    configuration.

    <location> : path of the directory where parfiles are found.
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



# Main class ============================================================================
class Interface:
    '''A class to hold global variables as attributes and give
    clear and concise structure to the main() function.'''

    @wait_for_ok("parsing input")
    def __init__(self, config_file,
                 nums: int = None, # or any int-returning iterable
                 output_dir: Path = Path.cwd(),
                 dust_bin_mode: str = "auto",
                 read_gas_density=False,
                 read_gas_velocity=False,
                 settling=False,
                 mcfost_verbose=False):

        self.warnings = []
        # input checking
        if not isinstance(config_file, (str, Path)):
            raise TypeError(config_file)
        if not isinstance(output_dir, (str, Path)):
            raise TypeError(output_dir)

        # attribute storage
        self._base_args = {
            'config_file': Path(config_file),
            'output_dir': Path(output_dir),
            'nums': nums,
            'dust_bin_mode': dust_bin_mode,
            'read_gas_density': read_gas_density,
            'read_gas_velocity': read_gas_velocity
        }

        self._dim = 2  # no support for 3D input yet
        self.mcfost_verbose = mcfost_verbose
        self.read_gas_velocity = read_gas_velocity
        self.use_settling = settling

        # parse configuration file
        self.config = f90nml.read(config_file)
        if nums is None:
            nums = self.config["amrvac_input"]["nums"]
        if isinstance(nums, int):
            self.nums = [nums]  # make it iterable
        else:
            self.nums = list(set(nums))  #make it iterable, filter out duplicates and sort them
        self.current_num = self.nums[0]

        hydro_data_dir = shell_path(self.config["amrvac_input"]["hydro_data_dir"])
        if not hydro_data_dir.is_absolute():
            options = self.config['amrvac_input']
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
            self.config['amrvac_input'].update({'hydro_data_dir': p.resolve()})
        self.sim_conf = read_amrvac_parfiles(
            parfiles=self.config['amrvac_input']['config'],
            location=self.config['amrvac_input']['hydro_data_dir']
        )

        self._µsizes = None
        self._input_data = None
        self._output_grid = None
        self._new_2D_arrays = None
        self._new_3D_arrays = None
        self._dust_binning_mode = None

        self._set_dust_binning_mode(dust_bin_mode, warning=False)
        if dust_bin_mode == "auto":
            self._autoset_dbm()
            assert self.dust_binning_mode != "auto"

        if not self.io.OUT.directory.exists():
            os.makedirs(self.io.OUT.directory)
            self.warnings.append(f"dir {self.io.OUT.directory} was created")

        if not self.config.get("units"):
            self.warnings.append(f"&units parameter list not found. Assuming {DEFAULT_UNITS}")
            self.config["units"] = f90nml.Namelist(DEFAULT_UNITS)
        else:
            for k, v in DEFAULT_UNITS.items():
                if not self.config["units"].get(k):
                    self.warnings.append(f"&units:{k} parameter not found. Assuming default {v}")
                    self.config["units"][k] = v

    @property
    def io(self) -> IOinfo:
        """Give up-to-date information on data location and naming (.i: input, .o: output)"""
        vtu_filename = "".join([self.sim_conf["filelist"]["base_filename"],
                                str(self.current_num).zfill(4),
                                ".vtu"])

        geomdefs = {"nr": 1, "nphi": 2}
        _input = DataInfo(
            directory=shell_path(self.config["amrvac_input"]["hydro_data_dir"]).resolve(),
            filename=vtu_filename,
            gridshape=GridShape(**{k: self.sim_conf["meshlist"][f"domain_nx{n}"]
                                   for k, n in geomdefs.items()})
        )

        trad_keys = {"nr": "n_rad", "nphi": "n_az", "nz": "nz"}
        _output = DataInfo(
            directory=Path(self._base_args["output_dir"]),
            filename=_input.filestem+".fits",
            gridshape=GridShape(**{k1: self.config["mcfost_output"][k2]
                                   for k1, k2 in trad_keys.items()})
        )
        return IOinfo(IN=_input, OUT=_output)

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
        elif self._bin_gas():
            self.warnings.append(
                f"read_gas_density asked but redundant in '{self.dust_binning_mode}' mode, ignored")
            rgd = False
        else:
            rgd = True
        return rgd

    def display_warnings(self):
        """A colorful way to print the warning list."""
        print(" WARNINGS:")
        print(RED+"\n".join([f" - {w}" for w in self.warnings]))
        if colorama is not None:
            print(colorama.Style.RESET_ALL, end="")


    # dust binning mode API
    # ================================================================
    @property
    def dust_binning_mode(self):
        """Define binning strategy
        - (gas-only)  : use only gas as a proxy for dust
        - (dust-only) : use only dust information
        - (mixed)     : use both, assuming gas traces the smallest grains
        """
        return self._dust_binning_mode

    def _autoset_dbm(self) -> None:
        """From dust_binning_mode=="auto" mode, select the correct one"""
        try:
            smallest_gs_µm = 1e4* min(np.array(self.sim_conf['usr_dust_list']['grain_size_cm']))
        except KeyError:
            self._set_dust_binning_mode("gas-only", reason="could not find grain sizes")
        else:
            if smallest_gs_µm > MINGRAINSIZE_µ:
                self._set_dust_binning_mode(
                    "mixed", reason=f"smallest size found > {MINGRAINSIZE_µ}µm"
                )

    def _set_dust_binning_mode(self, new_dbm: str, reason: str = None, warning=True):
        """Set value and add a warning."""
        if new_dbm not in {"dust-only", "gas-only", "mixed", "auto"}:
            raise KeyError(f'Unknown dust binning mode "{new_dbm}"')

        if warning:
            w = ["dust-binning mode was switched"]
            old = self._dust_binning_mode
            if old is not None:
                w.append(f'''from "{old}"''')
            w.append(f'''to "{new_dbm}"''')
            if reason is not None:
                w.append(f"\n   REASON: {reason}")
            self.warnings.append(" ".join(w))
        self._dust_binning_mode = new_dbm

    def _bin_dust(self) -> bool:
        """Should dust fluids be passed to mcfost ?"""
        return self.dust_binning_mode in {"dust-only", "mixed"}

    def _bin_gas(self) -> bool:
        """Should gas be passed to mcfost ?"""
        return self.dust_binning_mode in {'gas-only', 'mixed'}

    @property
    def grain_micron_sizes(self) -> np.ndarray:
        """Read grain sizes (assumed in [cm]), from AMRVAC parameters and
        convert to microns."""
        assert self.dust_binning_mode != "auto"
        if self._µsizes is None:
            µm_sizes = np.empty(0)
            if self._bin_dust():
                µm_sizes = 1e4 * np.array(
                    self.sim_conf["usr_dust_list"]["grain_size_cm"])
                assert min(µm_sizes) > 0.1 #in case this triggers, review this code
            # always associate a grain size to the gas bin
            µm_sizes = np.insert(µm_sizes, 0, MINGRAINSIZE_µ)
            self._µsizes = µm_sizes
        return self._µsizes

    @property
    def mcfost_conf_file(self) -> Path:
        """Locate output configuration file for mcfost"""
        return self.io.OUT.directory / "mcfost_conf.para"

    def load_input_data(self, n: int = None) -> None:
        '''Use vtkvacreader.VacDataSorter to load AMRVAC data'''
        if n is not None:
            assert n in self.nums
            #reinit properties
            self._input_data = None
            self._output_grid = None
            self._new_2D_arrays = None
            self._new_3D_arrays = None
            self.current_num = n
        self._input_data = VacDataSorter(
            file_name=str(self.io.IN.filepath),
            shape=(self.io.IN.gridshape.nr, self.io.IN.gridshape.nphi)
        )

    @property
    def input_data(self):
        '''Load input simulation data'''
        if self._input_data is None:
            self.load_input_data(self.current_num)
        return self._input_data

    @property
    def output_grid(self) -> dict:
        '''Store info on 3D output grid specifications
        as vectors "v", and (r-phi)grids "g"'''
        if self._output_grid is None:
            if not self.mcfost_conf_file.is_file():
                self.write_mcfost_conf_file()
            target_grid = get_mcfost_grid(self)
            self._output_grid = {
                'array': target_grid,
                # (nr, nphi) 2D grids
                'rg': target_grid[0, :, 0, :],
                'phig': target_grid[2, :, 0, :],
                # (nr, nz) 2D grid (z points do not depend on phi)
                'zg': target_grid[1, 0, :, :],
                # vectors (1D arrays)
                'rv': target_grid[0, 0, 0, :],
                'phiv': target_grid[2, :, 0, 0],
            }
        return self._output_grid

    @property
    def g2d_ratio(self):
        """Gas to dust ratio"""
        res = 0.01
        try:
            res = self.sim_conf["usr_dust_list"]["gas2dust_ratio"]
        except KeyError:
            self.warnings.append(f"could not find &usr_dust_list:gas2dust_ratio, assume {res}")
        return res

    def estimate_dust_mass(self) -> float:
        """Estimate the total dust mass in the grid, in solar masses"""
        # devnote : this assumes a linearly spaced grid
        dphi = 2*np.pi / self.io.IN.gridshape.nphi
        rvect = self.input_data.get_ticks(0)
        dr = rvect[1] - rvect[0]
        cell_surfaces = dphi/2 * ((rvect + dr/2)**2 - (rvect - dr/2)**2)

        if self.dust_binning_mode == "gas-only":
            keys = ["rho"]
        else:
            keys = [k for k, _ in self.input_data if "rhod" in k]
        mass = 0.0
        for key in keys:
            mass += np.sum([cell_surfaces * self.input_data[key][:, i]
                            for i in range(self.io.IN.gridshape.nphi)])
        if self.dust_binning_mode == "gas-only":
            mass /= self.g2d_ratio
        mass *= self.config["units"]["mass2solar"]
        return mass

    def _translate_amrvac_config(self) -> dict:
        parameters = {}

        # Zone
        mesh = self.sim_conf['meshlist']
        conv2au = self.config["units"]["distance2au"]
        parameters.update({
            'rin': mesh['xprobmin1']*conv2au,
            'rout': mesh['xprobmax1']*conv2au,
            'maps_size': 2*mesh['xprobmax1']*conv2au,
        })

        if self._bin_dust(): #devnote : using a private method outside of class...
            parameters.update({
                "gas_to_dust_ratio": self.g2d_ratio,
                "dust_mass": self.estimate_dust_mass()
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
            parameters.update({"star_mass": self.sim_conf["disk_list"]["central_mass"]})
        except KeyError:
            self.warnings.append("&disk_list not found. Assuming default values")
        return parameters

    def write_mcfost_conf_file(self) -> None:
        '''Customize defaults with user specifications'''
        custom = {}
        custom.update(self._translate_amrvac_config())
        unknown_args = self._scan_for_unknown_arguments()
        if unknown_args:
            raise KeyError(f'Unrecognized MCFOST argument(s): {unknown_args}')
        custom.update(self.config['mcfost_output'])

        write_mcfost_conf(
            output_file=self.mcfost_conf_file,
            custom=custom,
            verbose=self.mcfost_verbose
        )

    def _scan_for_unknown_arguments(self) -> list:
        """Get unrecognized arguments found in mcfost_output"""
        unknowns = []
        for arg in self.config["mcfost_output"].keys():
            if not arg.lower() in KNOWN_MCFOST_ARGS:
                unknowns.append(arg)
        return unknowns

    def write_output(self) -> None:
        """Write a .fits file suited for MCFOST input."""
        dust_bin_selector = {
            "gas-only": np.zeros(1, dtype="int64"),
            "dust-only": 1 + self.grain_micron_sizes[1:].argsort(),
            "mixed": self.grain_micron_sizes.argsort()
        }[self.dust_binning_mode]

        suppl_hdus = []
        assert (len(dust_bin_selector) > 1) == (self._bin_dust())
        if len(dust_bin_selector) > 1:
            # mcfost requires an HDU with grain sizes only if more than one population is present
            suppl_hdus.append(
                fits.ImageHDU(self.grain_micron_sizes[dust_bin_selector])
            )

        header = {'read_n_a': 0} # automatic normalization of size-bins from mcfost param file.
        if self.read_gas_density:
            #devnote: add try statement here ?
            header.update(dict(gas_to_dust=self.sim_conf["usr_dust_list"]["gas2dust_ratio"]))
            suppl_hdus.append(fits.ImageHDU(self.new_3D_arrays[0]))
            header.update(dict(read_gas_density=1))

        if self.read_gas_velocity:
            header.update(dict(read_gas_velocity=1))
            suppl_hdus.append(fits.ImageHDU(self.new_3D_gas_velocity))

        dust_densities_HDU = fits.PrimaryHDU(self.new_3D_arrays[dust_bin_selector])
        for k, v in header.items():
            # this is the canonical way to avoid HIERARCH-related warnings from astropy
            if len(k) > 8:
                k = f"HIERARCH {k}"
            dust_densities_HDU.header.append((k, v))

        with open(self.io.OUT.filepath, mode="wb") as fo:
            hdul = fits.HDUList(hdus=[dust_densities_HDU] + suppl_hdus)
            hdul.writeto(fo)

    @property
    def input_grid(self) -> dict:
        """Store physical coordinates (vectors) about the input grid specifications."""
        ig = {
            "rv": self.input_data.get_ticks("r") * self.config["units"]["distance2au"],
            "phiv": self.input_data.get_ticks("phi")
        }
        return ig

    def _interpolate2D(self, datakey: str) -> np.ndarray:
        """Transform a polar field from MPI-AMRVAC coords to mcfost coords"""
        interpolator = interp2d(
            self.input_grid["phiv"],
            self.input_grid["rv"],
            self.input_data[datakey], kind="cubic"
        )
        return interpolator(self.output_grid["phiv"], self.output_grid["rv"])

    def gen_2D_arrays(self) -> None:
        """Interpolate input data density fields from input coords to output coords"""
        density_keys = sorted(filter(lambda k: "rho" in k, self.input_data.fields.keys()))
        self._new_2D_arrays = np.array([self._interpolate2D(datakey=k) for k in density_keys])
        assert self._new_2D_arrays[0].shape == (self.io.OUT.gridshape.nr,
                                                self.io.OUT.gridshape.nphi)

    @property
    def aspect_ratio(self):
        """Dimensionless gas scale height implied by mcfost parameters"""
        mcfl = self.config["mcfost_output"]
        return mcfl["scale_height"] / mcfl["reference_radius"]

    def gen_3D_arrays(self) -> None:
        """Interpolate input data onto full 3D output grid"""
        oshape = self.io.OUT.gridshape
        nr, nphi, nz = oshape.nr, oshape.nphi, oshape.nz

        nbins = len(self.new_2D_arrays)
        self._new_3D_arrays = np.zeros((nbins, nphi, nz, nr))
        for ir, r in enumerate(self.output_grid["rv"]):
            z_vect = self.output_grid["zg"][nz:, ir].reshape(1, nz)
            gas_height = r * self.aspect_ratio
            for i_bin, grain_µsize in enumerate(self.grain_micron_sizes):
                surface_density = self.new_2D_arrays[i_bin, ir, :]
                H = gas_height
                if self.use_settling:
                    H *= (grain_µsize / MINGRAINSIZE_µ)**(-0.5)
                gaussian = np.exp(-z_vect**2/ (2*H**2)) / (np.sqrt(2*np.pi) * H)
                self._new_3D_arrays[i_bin, :, :, ir] = \
                    gaussian * surface_density.reshape(nphi, 1)

    @property
    def new_2D_arrays(self) -> list:
        '''Last minute generation is used if required'''
        if self._new_2D_arrays is None:
            self.gen_2D_arrays()
        return self._new_2D_arrays

    @property
    def new_3D_arrays(self) -> list:
        '''Last minute generation is used if required'''
        if self._new_3D_arrays is None:
            self.gen_3D_arrays()
        return self._new_3D_arrays

    @property
    def new_3D_gas_velocity(self) -> np.ndarray:
        """Derive the 3D velocity field for gas velocity, in km/s"""
        rho, mr, mphi = map(self._interpolate2D, ["rho", "m1", "m2"])
        vr, vphi = map(lambda x: x/rho, [mr, mphi])
        phig = self.output_grid["phig"].transpose()
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
        conv = self.config["units"]
        dimvel = conv["distance2au"]*units.au / (conv["time2yr"]*units.yr)
        vel2kms = dimvel.to(units.m / units.s).value
        velarr = np.stack([vx, vy, vz], axis=3) * vel2kms
        return velarr.transpose()



# Verbose version of Interface ==========================================================
class VerbatimInterface(Interface):
    """A more talkative Interface"""
    @wait_for_ok(f"loading input data")
    def load_input_data(self, n: int = None) -> None:
        super().load_input_data(n)

    @wait_for_ok('writing mcfost configuration file')
    def write_mcfost_conf_file(self) -> None:
        super().write_mcfost_conf_file()

    @wait_for_ok('interpolating to mcfost grid')
    def gen_2D_arrays(self):
        super().gen_2D_arrays()

    @wait_for_ok('converting 2D arrays to 3D')
    def gen_3D_arrays(self):
        super().gen_3D_arrays()

    @wait_for_ok('building the .fits file')
    def write_output(self) -> None:
        super().write_output()



# Main function =========================================================================
def main(config_file: str,
         nums: int = None, # or any in-returning interable
         output_dir: str = '.',
         dust_bin_mode: str = "auto",
         read_gas_density=False,
         read_gas_velocity=False,
         settling=False,
         verbose=False,
         mcfost_verbose=False):
    '''Try to transform a .vtu file into a .fits'''
    print(decorated_centered_message("start vac2fost"))
    InterfaceType = {True: VerbatimInterface, False: Interface}[verbose]
    itf = InterfaceType(config_file, nums=nums, output_dir=output_dir,
                        dust_bin_mode=dust_bin_mode,
                        read_gas_density=read_gas_density,
                        read_gas_velocity=read_gas_velocity,
                        settling=settling,
                        mcfost_verbose=mcfost_verbose)

    for i, n in enumerate(itf.nums):
        if verbose or i == 0:
            print()
        mess1 = f"current input number: {n}"
        mess2 = f"({i+1}/{len(itf.nums)})"
        print((BOLD+" "*(get_prompt_size()-len(mess1)-len(mess2))).join([mess1, mess2]))
        print("-"*get_prompt_size())
        try:
            itf.load_input_data(n)
        except FileNotFoundError as err:
            filepath = Path(str(err)).relative_to(Path.cwd())
            itf.warnings.append(f"file not found: {filepath}")
            continue
        itf.write_mcfost_conf_file()
        itf.gen_2D_arrays()
        itf.gen_3D_arrays()
        itf.write_output()

        try:
            filepath = itf.io.OUT.filepath.relative_to(Path.cwd())
        except ValueError:
            filepath = itf.io.OUT.filepath
        print(CYAN + f" >>> wrote {filepath}")

    if itf.warnings:
        print()
        itf.display_warnings()

    print(decorated_centered_message("end vac2fost"))

    # return the Interface object for inspection (tests)
    return itf



# Script part ===========================================================================
if __name__ == '__main__':
    # Parse the script arguments
    parser = ArgumentParser(description='Parse arguments for main app')
    parser.add_argument(
        dest='configuration', type=str,
        nargs='?',
        default=None,
        help='configuration file (namelist) for this script'
    )
    parser.add_argument(
        "-n", "--nums", dest="nums", type=int,
        required=False,
        default=None,
        nargs="*",
        help="output number(s) of the target .vtu VAC output file to be converted"
    )
    parser.add_argument(
        '-o', '--output', dest='output', type=str,
        required=False,
        default='.',
        help='select output directory for generated files'
    )
    parser.add_argument(
        '-dbm', '--dustbinmode', dest='dbm', type=str,
        required=False,
        default="auto",
        help="prefered bin selection mode [dust-only, gas-only, mixed, auto]"
    )
    parser.add_argument(
        "--read_gas_density",
        action="store_true",
        help="pass gas density to mcfost"
    )
    parser.add_argument(
        "--read_gas_velocity",
        action="store_true",
        help="pass gas velocity to mcfost (keplerian velocity is assumed otherwise)"
    )
    parser.add_argument(
        "--settling",
        action="store_true",
        help="differentiate scale height according to grain-size"
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='activate verbose mode'
    )
    parser.add_argument(
        '--mcfost_verbose',
        action='store_true',
        help='do not silence mcfost'
    )
    parser.add_argument(
        '--genconf', action='store_true',
        help="print a default configuration file for vac2fost"
    )
    parser.add_argument(
        '--cprofile',
        action='store_true',
        help='activate code profiling'
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="display this code's version"
    )

    cargs = parser.parse_args()

    if cargs.version:
        print(__version__)
        sys.exit(0)

    if cargs.genconf:
        print(generate_conf_template())
        print(f"%% automatically generated with vac2fost {__version__}\n")
        sys.exit(0)
    elif len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if cargs.cprofile:
        import cProfile
        import pstats
        import io
        pr = cProfile.Profile()
        pr.enable()
    # -------------------------------------------
    main(
        config_file=cargs.configuration,
        nums=cargs.nums,
        output_dir=cargs.output.strip(),
        dust_bin_mode=cargs.dbm,
        read_gas_density=cargs.read_gas_density,
        read_gas_velocity=cargs.read_gas_velocity,
        settling=cargs.settling,
        verbose=cargs.verbose,
        mcfost_verbose=cargs.mcfost_verbose
    )
    # -------------------------------------------
    if cargs.cprofile:
        pr.disable()
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
        ps.print_stats()
        print(s.getvalue())
