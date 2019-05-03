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
__version__ = "2.3.3"
min_mcfost_version = "3.0.35"  # minimal requirement
#rec_mcfost_version = "3.0.35"  # recommendation



# Imports
# =======================================================================================
# stdlib
import os
import sys
import shutil
from dataclasses import dataclass
from pathlib import Path
from subprocess import run, CalledProcessError
from argparse import ArgumentParser
from collections import OrderedDict as od
from socket import gethostname
from uuid import uuid4 as uuid

# non standard externals
import numpy as np
from astropy import units
from astropy.io import fits
from scipy.interpolate import interp2d
import f90nml
try:
    import colorama
    colorama.init(autoreset=True)
    BOLD = colorama.Style.BRIGHT
    RED = BOLD + colorama.Fore.RED
    CYAN = colorama.Fore.CYAN
except ImportError:
    colorama = None
    BOLD = RED = CYAN = ""

# private externals
from vtk_vacreader import VacDataSorter



# mcfost detection ======================================================================
if shutil.which("mcfost") is None:
    raise OSError(RED+"could not find mcfost. Please install mcfost before using vac2fost")
out = run("yes | mcfost -version", shell=True, capture_output=True).stdout #binary
out = "".join(map(chr, out))

DETECTED_MCFOST_VERSION = out.split("\n")[0].split()[-1]
del out
if DETECTED_MCFOST_VERSION < min_mcfost_version:
    raise OSError(f"mcfost version must be >= {min_mcfost_version}")



# Globals ===============================================================================
MINGRAINSIZE_µ = 0.1
DEFAULT_UNITS = dict(distance2au=1.0, time2yr=1.0, mass2solar=1.0)



# Defintions ============================================================================
@dataclass
class GridShape:
    """Describe number of cells in cylindrical coordinates in a grid"""
    nr: int
    nphi: int
    nz: int = 1

@dataclass
class DataInfo:
    """Hold basic info about input or output data location and shape"""
    directory: Path
    filename: str
    gridshape: GridShape

    @property
    def filepath(self) -> Path:
        """full path"""
        return self.directory / self.filename

    @property
    def filestem(self) -> str:
        """filename without an extension (suffix)"""
        return str(Path(self.filename).stem)

@dataclass
class IOinfo:
    """Hold input and output data information"""
    IN: DataInfo
    OUT: DataInfo

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
        ref_radius=100.0,  # [a.u.]
        scale_height=1.0,  # [a.u.], at defined at ref_radius
    )
    sublists = {
        "amrvac_input": amrvac_list,
        "units": DEFAULT_UNITS,
        "mcfost_output": mcfost_list
    }
    template = f90nml.Namelist({k: f90nml.Namelist(v) for k, v in sublists.items()})
    return template

# decorators
def get_prompt_size():
    """size of command line interface messages sized to window. Caps at 80."""
    cols, _ = shutil.get_terminal_size()
    return min(cols, 80)

def parameterized(dec):
    """meta decorator, allow definition of decorators with parameters
    source: https://stackoverflow.com/questions/5929107/decorators-with-parameters"""
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)
        return repl
    return layer

@parameterized
def wait_for_ok(func, mess, lenght=get_prompt_size()-7):
    """decorator, sandwich the function execution with '<mess>  ...' & 'ok'"""
    def modfunc(*args, **kwargs):
        print(mess.ljust(lenght), end="... ", flush=True)
        result = func(*args, **kwargs)
        print("ok")
        return result
    return modfunc

def shell_path(pin: str) -> Path:
    """Transform <pin> to a Path, expanding included env variables."""
    return Path(os.path.expandvars(str(pin)))

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

class MCFOSTUtils:
    """Utility functions to call MCFOST in vac2fost.main()
    to define the output grid."""

    blocks_descriptors = od(
        # This nested orderded dictionnary describes a default parafile for mcfost
        #
        # parameter names should match mcfost's documentation
        # http://ipag-old.osug.fr/~pintec/mcfost/docs/html/parameter_file.html
        #
        # this is still WIP (remaining naming discrepencies)
        #
        # notes for a future pull-request on mcfost documentation itself:
        # *   indicates a change in name for various reasons...
        # ?   indicates a missing documentation line (or a deprecated parameter)
        # $   indicates stuff I'll have to go over again, either because it breaks regression here, or because I'll need to change to api altogether

        # fixplan : DONE 1) PR to fix <> and % in different commits (and %% ?)
        #                2) ask Christophe about "?" and go over *
        #                3) deal with *
        #                4) deal with $

        [
            ("Photons", (
                od([("nbr_photons_eq_temp", "1.28e5")]),
                od([("nbr_photons_lambda", "1.28e3")]),
                od([("nbr_photons_image", "1.28e5")])
            )),
            ("Wavelengths", (
                od([("n_lambda", 50),
                    ("lambda_min", "0.1"),
                    ("lambda_max", "3e3")]),
                od([("ltemp", True),
                    ("lsed", True),
                    ("use_default_wavelength_grid", True)]),
                od([("wavelength_file", "wavelengths.dat")]),
                od([("separate_contributions", False),
                    ("output_stokes_parameters", False)])
            )),
            ("Grid", (
                od([("grid_type", 1)]),
                od([("n_rad", 100),
                    ("nz", 10),
                    ("n_az", 100),
                    ("n_rad_in", 30)])
            )),
            ("Images", (
                od([("grid_nx", 501),
                    ("grid_ny", 501),
                    ("map_size", 400)]),
                od([("RT_imin", 0),
                    ("RT_imax", 0),
                    ("RT_n_incl", 1),
                    ("RT_centered", False)]),
                od([("RT_az_min", 0),
                    ("RT_az_max", 240),
                    ("RT_n_az", 1)]),
                od([("distance", 140)]),
                od([("disk_PA", 0)])
            )),
            ("Scattering Method", (
                od([("scattering_method", 0)]),
                od([("Mie_hg", 1)]) # *
            )),
            ("Symmetries", (
                od([("image_symmetry", True)]),
                od([("central_symmetry", True)]),
                od([("plane_symmetry", True)]),
            )),
            ("Disk physics", (
                od([("dust_settling", 0),
                    ("exp_strat", 0.5),
                    ("a_srat", 1.0)]),
                od([("dust_radial_migration", False)]),
                od([("sublimate_dust", False)]), # TODO: check order !!
                od([("hydrostatic_equilibrium", False)]),
                od([("viscous_heating", False),
                    ("viscosity", "1e-3")]),
            )),
            ("Number of Zones", (
                od([("n", 1)]),
            )),
            ("Density structure", (
                od([("zone_type", 1)]),
                od([("disk_dust_mass", "1e-3"),
                    ("gas_to_dust_ratio", 100)]),
                od([("scale_height", 5.0), # $
                    ("ref_radius", 100.0), # $
                    ("vertical_profile_exponent", 2)]),
                # todo : this part is oddly hard to rewrite... (breaks regression)
                od([('rin', 10),   #$
                    ('edge', 0),
                    ('rout', 200), #$
                    ('Rc', 100)]),
                # ^^^^^^^^^^^
                od([("flaring_exp", 1.0)]), # *
                od([("density_exp", -0.5), # *
                    ("gamma_exp", 0.0)]) # *
            )),
            ("Grain properties", (
                od([("n_species", 1)]),
                od([("Grain_type", "Mie"),
                    ("n_components", 1),
                    ("mixing_rule", 2),
                    ("porosity", 0.),
                    ("mass_fraction", 0.75),
                    ("DHS_Vmax", 0.9)]),
                od([("optical_indices_file", "Draine_Si_sUV.dat"),
                    ("volume_fraction", 1.0)]),
                od([("heating_method", 1)]),
                od([("amin", MINGRAINSIZE_µ),
                    ("amax", 1000),
                    ("aexp", 3.5),
                    ("n_grains", 100)])
            )),
            ("Molecular RT settings", (
                od([("lpop", True),
                    ("lpop_accurate", True),
                    ("LTE", True),
                    ("profile_width", 15.0)]),
                od([("v_turb", 0.0)]), # ?
                od([("nmol", 1)]), # ?
                od([("molecular_data_file", "13co.dat"), # *
                    ("level_max", 6)]),
                od([("vmax", 1.0),
                    ("n_speed", 20)]),
                od([("cst_abundance", True),
                    ("abund", "1e-6"), # ?
                    ("abund_file", "abundance.fits.gz")]), # ?
                od([("ray_tracing", True), # *
                    ("n_lines", 3)]),
                od([("transition_num_1", 1),
                    ("transition_num_2", 2),
                    ("transition_num_3", 3)])
            )),
            ("Star properties", (
                od([("n_stars", 1)]),
                od([("Teff", 4000.0),
                    ("Rstar", 2.0),
                    ("Mstar", 1.0), # *
                    ("x", 0.),
                    ("y", 0.),
                    ("z", 0),
                    ("is_blackbody", True)]),
                od([("star_rad_file", "lte4000-3.5.NextGen.fits.gz")]), # ?
                od([("fUV", 0.0), ("slope_fUV", 2.2)]),
            ))
        ])

    known_args = []
    for descriptor in blocks_descriptors.items():
        for di in descriptor[1]:
            known_args += [k.lower() for k in di.keys()]

    def write_mcfost_conf(output_file: Path, custom: dict = None, verbose=False):
        """Write a configuration file for mcfost using values from <custom>,
        and falling back to defaults found in block_descriptor defined above
        """
        if custom is None:
            custom = {}
        if output_file.exists() and verbose:
            print(f'Warning: {output_file} already exists, and will be overwritten.')
        with open(output_file, mode="wt") as fi:
            fi.write(".".join(min_mcfost_version.split(".")[:2]).ljust(10) +
                     "mcfost minimal version prescribed by vac2fost\n\n")
            for block, lines in __class__.blocks_descriptors.items():
                fi.write(f'# {block}\n')
                for line in lines:
                    parameters = []
                    for param, default in line.items():
                        if param in custom:
                            val = custom[param]
                        else:
                            val = default
                        parameters.append(str(val))
                    fi.write('  ' + '  '.join(parameters).ljust(36)
                             + '  ' + ', '.join(line.keys()))
                    fi.write('\n')
                fi.write('\n')
            fi.write("\n\n")
            fi.write(f"%% automatically generated with vac2fost {__version__}\n")
            fi.write(f"%% via mcfost {DETECTED_MCFOST_VERSION}\n")
            fi.write(f"%% run by {os.environ['USER']} on {gethostname()}\n")
        if verbose:
            print(f'wrote {output_file}')

    def get_mcfost_grid(itf) -> np.ndarray:
        '''Pre-run MCFOST with -disk_struct flag to get the exact grid used.'''
        mcfost_conf_file = itf.mcfost_conf_file
        output_dir = itf.io.OUT.directory

        output_dir = Path(output_dir).resolve()
        mcfost_conf_path = Path(mcfost_conf_file)
        if not output_dir.exists():
            os.makedirs(output_dir)

        grid_file_name = output_dir / 'mcfost_grid.fits.gz'

        if itf.current_num == itf.nums[0]:
            assert mcfost_conf_path.exists()
            # generate a grid data file with mcfost itself and extract it
            tmp_mcfost_dir = output_dir / f"TMP_VAC2FOST_MCFOST_GRID_{uuid()}"
            os.makedirs(tmp_mcfost_dir)
            try:
                shutil.copyfile(mcfost_conf_path.resolve(),
                                tmp_mcfost_dir/mcfost_conf_path.name)
            except shutil.SameFileError:
                pass

            pile = Path.cwd()
            os.chdir(tmp_mcfost_dir)
            try:
                run(["mcfost", "mcfost_conf.para", "-disk_struct"], check=True,
                    capture_output=(not itf.mcfost_verbose))

                shutil.move("data_disk/grid.fits.gz", grid_file_name)
            except CalledProcessError as exc:
                errtip = f"\nError in MCFOST, exited with exitcode {exc.returncode}"
                if exc.returncode == 174:
                    errtip += (
                        "\nThis is probably a memory issue. "
                        "Try reducing the target resolution or,"
                        " alternatively, give more cpu memory to this task."
                    )
                raise RuntimeError(errtip) from exc
            finally:
                os.chdir(pile)
                shutil.rmtree(tmp_mcfost_dir)
        with fits.open(grid_file_name, mode='readonly') as fi:
            target_grid = fi[0].data
        return target_grid



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
        '''Read grain sizes (assumed in [cm]), from AMRVAC parameters and
        convert to microns.'''
        assert self.dust_binning_mode != "auto"
        if self._µsizes is None:
            µm_sizes = np.empty(0)
            if self._bin_dust():
                cm_sizes = np.array(
                    self.sim_conf['usr_dust_list']['grain_size_cm'])
                µm_sizes = 1e4 * cm_sizes
            if self._bin_gas():
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
            target_grid = MCFOSTUtils.get_mcfost_grid(self)
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

        MCFOSTUtils.write_mcfost_conf(
            output_file=self.mcfost_conf_file,
            custom=custom,
            verbose=self.mcfost_verbose
        )

    def _scan_for_unknown_arguments(self) -> list:
        """Get unrecognized arguments found in mcfost_output"""
        unknowns = []
        for arg in self.config["mcfost_output"].keys():
            if not arg.lower() in MCFOSTUtils.known_args:
                unknowns.append(arg)
        return unknowns

    def write_output(self) -> None:
        """Write a .fits file suited for MCFOST input."""
        dust_bin_selector = {
            "gas-only": np.zeros(1, dtype="int64"),
            "dust-only": 1 + self.grain_micron_sizes.argsort(),
            "mixed": self.grain_micron_sizes.argsort()
        }[self.dust_binning_mode]

        suppl_hdus = []
        if len(dust_bin_selector) > 1:
            # mcfost requires an HDU with grain sizes only if more than one population is present
            suppl_hdus.append(
                fits.ImageHDU(self.grain_micron_sizes[self.grain_micron_sizes.argsort()])
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
        """Dimensionless ratio implied by mcfost parameters"""
        mcfl = self.config['mcfost_output']
        return mcfl['scale_height'] / mcfl['ref_radius']

    def gen_3D_arrays(self) -> None:
        """Interpolate input data onto full 3D output grid"""
        oshape = self.io.OUT.gridshape
        nr, nphi, nz = oshape.nr, oshape.nphi, oshape.nz

        nbins = len(self.new_2D_arrays)
        self._new_3D_arrays = np.zeros((nbins, nphi, nz, nr))
        for ir, r in enumerate(self.output_grid["rv"]):
            z_vect = self.output_grid["zg"][nz:, ir].reshape(1, nz)
            local_height = r * self.aspect_ratio
            gaussian = np.exp(-z_vect**2/ (2*local_height**2)) / (np.sqrt(2*np.pi) * local_height)
            for i_bin, surface_density in enumerate(self.new_2D_arrays[:, ir, :]):
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



# Main function and associated prompt utils =============================================
def decorated_centered_message(mess: str, dec: str = "=") -> str:
    """Return a decorated version of <mess>"""
    ndecor = int((get_prompt_size() - (len(mess)+2)) / 2)
    return BOLD + " ".join([dec*ndecor, mess, dec*ndecor])

def main(config_file: str,
         nums: int = None, # or any in-returning interable
         output_dir: str = '.',
         dust_bin_mode: str = "auto",
         read_gas_density=False,
         read_gas_velocity=False,
         verbose=False,
         mcfost_verbose=False):
    '''Try to transform a .vtu file into a .fits'''
    print(decorated_centered_message("start vac2fost"))
    InterfaceType = {True: VerbatimInterface, False: Interface}[verbose]
    itf = InterfaceType(config_file, nums=nums, output_dir=output_dir,
                        dust_bin_mode=dust_bin_mode,
                        read_gas_density=read_gas_density,
                        read_gas_velocity=read_gas_velocity,
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
