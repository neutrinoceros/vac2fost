#!/usr/bin/env python3
'''A conversion tool from AMRVAC output (.vtu) to MCFOST input (.fits)

Run `vac2fost.py --help` for documentation on command line usage

Steps
  a) load AMRVAC data with vtk_vacreader.VacDataSorter(), sort it as 2D arrays
  b) dry-run MCFOST to get the exact target grid
  c) interpolate data to the target grid
  d) convert to 3D (gaussian redistribution of density)
  e) collect, sort and write output data to a fits file

Known limitations
  1) amr is not supported by reader
  2) portability is not guaranted
  3) interpolation does not account for the curvature of polar cells
  4) only r-phi input grids are currently supported
'''
__version__ = "2.3.0"
mcfost_major_version = "3.0"
mcfost_minor_version = "35"



# Imports
# =======================================================================================
# stdlib
import os
import sys
import shutil
from dataclasses import dataclass
from warnings import warn
from pathlib import Path
from subprocess import run, CalledProcessError
from argparse import ArgumentParser
from collections import OrderedDict as od
from socket import gethostname
from uuid import uuid4 as uuid

# non standard externals
import numpy as np
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
try:
    run(["which", "mcfost"], check=True, capture_output=True)
except CalledProcessError:
    print(RED+"Critical: could not find mcfost. Please install mcfost before using vac2fost")

bout = run("yes | mcfost -version", shell=True, capture_output=True).stdout
out = "".join(map(chr, bout))
version_tag = out.split("\n")[0].split()[-1]

verx, very, verz = map(int, version_tag.split("."))

if float(f"{verx}.{very}") < float(mcfost_major_version):
    raise EnvironmentError("mcfost version must be >= {mcfost_major_version}")

EXPECTED_ZSHAPE_INCREMENT = 0
if f"{verx}.{very}" == "3.0":
    if verz < 32:
        warn("vac2fost has not been tested for mcfost < 3.0.32")
    if verz < 35:
        EXPECTED_ZSHAPE_INCREMENT = 1

DETECTED_MCFOST_VERSION = verx, very, verz
del bout, out, version_tag, verx, very, verz



# Globals ===============================================================================
MINGRAINSIZE_µ = 0.1
S2YR = 1/(365*24*3600)
AU2KM = 149597870.700



# Defintions ============================================================================
@dataclass
class DataInfo:
    """Hold basic info about input or output data location and shape"""
    directory: Path
    filename: str
    shape: tuple

    @property
    def filepath(self):
        """full path"""
        return self.directory / self.filename

    @property
    def filestem(self):
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

    amrvac_unit_list = dict(distance2au=1.0, time2yr=1.0)

    mcfost_list = dict(
        nr=128, nr_in=4, nphi=128, nz=10,
        # aspect ratio is implied by those parameters
        flaring_index=1.125,
        ref_radius=100.0,  # [a.u.]
        scale_height=1.0,  # [a.u.], at defined at ref_radius
    )
    sublists = {
        "amrvac_input": amrvac_list,
        "units": amrvac_unit_list,
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

    confs = [f90nml.read(pathloc / pf) for pf in pfs]
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
        # name every mcfost parameter (by order of appearence)
        # and give default values
        [
            ('Photons', (
                od([('nphot_temp', '1.28e5')]),
                od([('nphot_sed', '1.28e3')]),
                od([('nphot_img', '1.28e5')])
            )),
            ('Wavelengths', (
                od([('n_lambda', 50),
                    ('lambda_min', 0.1),
                    ('lambda_max', 3e3)]),
                od([('compute_temp', True),
                    ('compute_sed', True),
                    ('use_default_wl', True)]),
                od([('wavelength_file',
                     'wavelengths.dat')]),
                od([('separation', False),
                    ('stokes_parameters', False)])
            )),
            ('Grid', (
                od([('geometry', '1')]),
                od([('nr', 100),
                    ('nz', 10),
                    ('nphi', 100),
                    ('nr_in', 30)])
            )),
            ('Maps', (
                od([('nx', 501),
                    ('ny', 501),
                    ('maps_size', 400)]),
                od([('imin', 0),
                    ('imax', 0),
                    ('n_incl', 1),
                    ('centered', False)]),
                od([('az_min', 0),
                    ('az_max', 240),
                    ('n_az_angles', 1)]),
                od([('distance_pc', 140)]),
                od([('disk_position_angle', 0)])
            )),
            ('Scattering', (
                od([('scattering_method', '0')]),
                od([('theory', 1)])
            )),
            ('Symmetries', (
                od([('sym_image', False)]),
                od([('sym_central', False)]),
                od([('sym_axial', False)]),
            )),
            ('Disk physics', (
                od([('dust_settling', 3),
                    ('exp_strat', 0.5),
                    ('a_srat', 1.0)]),
                od([('dust_radial_migration', False)]),
                od([('sublimate_dust', False)]),
                od([('hydrostatic_eq', False)]),
                od([('viscous_heating', False),
                    ('alpha_viscosity', '1e-3')]),
            )),
            ('Number of zones', (
                od([('n_zones', '1')]),
            )),
            ('Zone', (
                od([('zone_type', 1)]),
                od([('dust_mass', '1e-3'),
                    ('gas_to_dust_ratio', 100)]),
                od([('scale_height', 5.0),
                    ('ref_radius', 100.0),
                    # only relevant in geometry "4" (debris disk) according to doc
                    ('profile_exp', 2)]),
                od([('rin', 10),
                    ('edge', 0),
                    ('rout', 200),
                    ('rc', 100)]),
                od([('flaring_index', 1.0)]),
                od([('density_exp', -0.5),
                    ('gamma_exp', 0.0)])
            )),
            ('Grains', (
                od([('n_species', 1)]),
                od([('grain_type', 'Mie'),
                    ('n_components', 1),
                    ('mixing_rule', 2),
                    ('porosity', 0.),
                    ('mass_fraction', 0.75),
                    ('vmax_dhs', 0.9)]),
                od([('optical_indices_file', 'Draine_Si_sUV.dat'),
                    ('volume_fraction', 1.0)]),
                od([('heating_method', 1)]),
                od([('sp_min', MINGRAINSIZE_µ),
                    ('sp_max', 1000),
                    ('sexp', 3.5),
                    ('n_grains', 100)])
            )),
            ('Molecular RT', (
                od([('lpop', True),
                    ('laccurate_pop', True),
                    ('LTE', True),
                    ('profile_width', 15.)]),
                od([('v_turb', 0.2)]),
                od([('nmol', 1)]),
                od([("mol_data_file", "13co.dat"),
                    ('level_max', 6)]),
                od([('vmax', 1.0),
                    ('n_speed', 20)]),
                od([('cst_mol_abund', True),
                    ('abund', '1e-6'),
                    ('abund_file', 'abundance.fits.gz')]),
                od([('ray_tracing', True),
                    ('n_lines_rt', 3)]),
                od([('transition_num_1', 1),
                    ('transition_num_2', 2),
                    ('transition_num_3', 3)])
            )),
            ('Star', (
                od([('n_stars', 1)]),
                od([('star_temp', 4000.0),
                    ('star_radius', 2.0),
                    ('star_mass', 1.0),
                    ('star_x', 0.),
                    ('star_y', 0.),
                    ('star_z', 0),
                    ('star_is_bb', True)]),
                od([('star_rad_file', 'lte4000-3.5.NextGen.fits.gz')]),
                od([('fUV', 0.0), ('slope_fUV', 2.2)]),
            ))
        ])

    known_args = []
    for descriptor in blocks_descriptors.items():
        for di in descriptor[1]:
            known_args += list(di.keys())

    def write_mcfost_conf(output_file: Path, custom: dict = None, verbose=False):
        """Write a configuration file for mcfost using values from <custom>,
        and falling back to defaults found in block_descriptor defined above
        """
        if custom is None:
            custom = {}
        if output_file.exists() and verbose:
            print(f'Warning: {output_file} already exists, and will be overwritten.')
        with open(output_file, mode="wt") as fi:
            fi.write(mcfost_major_version.ljust(10) +
                     f"mcfost minimal version. Recommended minor {mcfost_minor_version}\n\n")
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
            fi.write(f"%% run by {os.environ['USER']} on {gethostname()}\n")
        if verbose:
            print(f'wrote {output_file}')

    def translate_amrvac_config(itf) -> dict:
        # itf must be of type Interface (can't be parsed properly before python 3.7)
        # devnote : this should be refactored as part of the Interface class
        '''pass amrvac parameters to mcfost'''
        parameters = {}

        # Zone
        mesh = itf.sim_conf['meshlist']
        conv2au = itf.config["units"]["distance2au"]
        parameters.update({
            'rin': mesh['xprobmin1']*conv2au,
            'rout': mesh['xprobmax1']*conv2au,
            'maps_size': 2*mesh['xprobmax1']*conv2au,
        })

        if itf._bin_dust(): #devnote : using a private method outside of class...
            dl2 = itf.sim_conf['usr_dust_list']
            parameters.update({
                'gas_to_dust_ratio': dl2['gas2dust_ratio'],
                # 'dust_mass': ... # MISSING FEATURE (ISSUE 18)
            })
            # Grains
            sizes_µm = itf.grain_micron_sizes
            parameters.update({
                # min/max grain sizes in microns
                'sp_min': min(1e-1, min(sizes_µm)),
                'sp_max': max(1e3, max(sizes_µm)),
            })
        return parameters

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
                print(errtip)
                raise
            finally:
                os.chdir(pile)
                shutil.rmtree(tmp_mcfost_dir)
        with fits.open(grid_file_name, mode='readonly') as fi:
            target_grid = fi[0].data
        return target_grid


def get_dust_mass(data: VacDataSorter) -> float:
    '''estimate the total dust mass in the grid in code units
    (solar mass = 1) is assumed by the present script and MCFOST
    '''
    # devnote : assume a linearly spaced grid
    dphi = 2*np.pi / data.shape[1]
    rvect = data.get_ticks(0)
    dr = rvect[1] - rvect[0]
    cell_surfaces = dphi/2 * ((rvect + dr/2)**2 - (rvect - dr/2)**2)

    mass = 0.0
    for _, field in filter(lambda item: 'rhod' in item[0], data):
        mass += np.sum([cell_surfaces * field[:, i]
                        for i in range(field.shape[1])])
    return mass


# Main class ============================================================================
class Interface:
    '''A class to hold global variables as attributes and give
    clear and concise structure to the main() function.'''

    @wait_for_ok("parsing input")
    def __init__(self, config_file,
                 nums: int = None, # or any int-returning iterable
                 output_dir: Path = Path('.'),
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

            elif not any(found):
                raise FileNotFoundError(hydro_data_dir/options['config'][0])
            else:
                p = (p1, p2)[found.index(True)]
            self.config['amrvac_input'].update({'hydro_data_dir': p.resolve()})
        self.sim_conf = read_amrvac_parfiles(
            parfiles=self.config['amrvac_input']['config'],
            location=self.config['amrvac_input']['hydro_data_dir']
        )

        self._dust_binning_mode = dust_bin_mode
        if dust_bin_mode == "auto":
            self._autoset_dbm()
            assert self.dust_binning_mode != "auto"

        self._µsizes = None

        self._input_data = None
        self._output_grid = None
        self._new_2D_arrays = None
        self._new_3D_arrays = None

        if not self.io.OUT.directory.exists():
            os.makedirs(self.io.OUT.directory)
            self.warnings.append(f"dir {self.io.OUT.directory} was created")

        # optional definition of the distance unit
        default_units = dict(distance2au=1.0, time2yr=1.0)
        if not self.config.get("units"):
            self.warnings.append(f"&units parameter list not found. Assuming {default_units}")
            self.config["units"] = f90nml.Namelist(default_units)
        else:
            for k, v in default_units.items():
                if not self.config["units"].get(k):
                    self.warnings.append(f"&units:{k} parameter not found. Assuming default {v}")
                    self.config["units"][k] = v

    @property
    def io(self) -> IOinfo:
        """Give up-to-date information on data location and naming (.i: input, .o: output)"""
        vtu_filename = "".join([self.sim_conf["filelist"]["base_filename"],
                                str(self.current_num).zfill(4),
                                ".vtu"])

        _input = DataInfo(
            directory=shell_path(self.config["amrvac_input"]["hydro_data_dir"]).resolve(),
            filename=vtu_filename,
            shape=tuple(
                [self.sim_conf["meshlist"][f"domain_nx{n}"] for n in range(1, self._dim+1)])
        )

        outshape = self.config["mcfost_output"]
        _output = DataInfo(
            directory=Path(self._base_args["output_dir"]),
            filename=_input.filestem+".fits",
            shape=(outshape["nr"], outshape["nz"], outshape["nphi"])
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

    def _set_dust_binning_mode(self, new_dbm: str, reason: str = None):
        """Set value and add a warning."""
        if new_dbm not in {"dust-only", "gas-only", "mixed", "auto"}:
            raise KeyError(f'Unknown dust binning mode "{new_dbm}"')

        old = self._dust_binning_mode
        w = f'dust-binning mode was switched from "{old}" to "{new_dbm}"'
        if reason is not None:
            w += f"\n   REASON: {reason}"
        self.warnings.append(w)
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
            shape=self.io.IN.shape
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

    def write_mcfost_conf_file(self) -> None:
        '''Customize defaults with user specifications'''
        custom = {}
        custom.update(MCFOSTUtils.translate_amrvac_config(self))
        unknown_args = self.scan_for_unknown_arguments()
        if unknown_args:
            raise KeyError(f'Unrecognized MCFOST argument(s): {unknown_args}')
        custom.update(self.config['mcfost_output'])

        if self._bin_dust():
            custom.update({'dust_mass': get_dust_mass(self.input_data)})
        MCFOSTUtils.write_mcfost_conf(
            output_file=self.mcfost_conf_file,
            custom=custom,
            verbose=self.mcfost_verbose
        )

    def scan_for_unknown_arguments(self) -> list:
        '''Get unrecognized arguments found in mcfost_output'''
        unknowns = []
        for arg in self.config['mcfost_output'].keys():
            if not arg.lower() in MCFOSTUtils.known_args:
                unknowns.append(arg)
        return unknowns

    def write_output(self) -> None:
        """Write a .fits file suited for MCFOST input."""
        dust_bin_selector = {
            "gas-only": 0,
            "dust-only": 1 + self.grain_micron_sizes.argsort(),
            "mixed": self.grain_micron_sizes.argsort()
        }[self.dust_binning_mode]

        additional_hdus = [
            fits.ImageHDU(self.grain_micron_sizes[self.grain_micron_sizes.argsort()])
        ]
        header = {'read_n_a': 0} # automatic normalization of size-bins from mcfost param file.
        if self.read_gas_density:
            #devnote: add try statement here ?
            header.update(dict(gas_to_dust=self.sim_conf["usr_dust_list"]["gas2dust_ratio"]))
            additional_hdus.append(fits.ImageHDU(self._new_3D_arrays[0]))
            header.update(dict(read_gas_density=1))

        if self.read_gas_velocity:
            header.update(dict(read_gas_velocity=1))
            rho, mr, mphi = map(self._interpolate2D, ["rho", "m1", "m2"])
            vr, vphi = map(lambda x: x/rho, [mr, mphi])
            phig = self.output_grid["phig"].transpose()
            vx = vr * np.cos(phig) + vphi * np.sin(phig)
            vy = vr * np.sin(phig) + vphi * np.cos(phig)

            # transform to 3D
            nz = self.io.OUT.shape[1]
            vx, vy = map(lambda a: np.stack([a]*nz, axis=1), [vx, vy])
            vz = np.zeros(vx.shape)

            # unit conversion
            units = self.config["units"]
            vel2km_per_s = units["distance2au"]*AU2KM / (units["time2yr"]*S2YR)
            vx *= vel2km_per_s
            vy *= vel2km_per_s

            # append
            for v in (vx, vy, vz):
                np.testing.assert_array_equal(v.shape, self.io.OUT.shape)

            additional_hdus.append(fits.ImageHDU(np.stack([vx, vy, vz], axis=3).T))

        dust_densities_HDU = fits.PrimaryHDU(self.new_3D_arrays[dust_bin_selector])
        for k, v in header.items():
            # this is the canonical way to avoid HIERARCH-related warnings from astropy
            if len(k) > 8:
                k = f"HIERARCH {k}"
            dust_densities_HDU.header.append((k, v))

        with open(self.io.OUT.filepath, mode="wb") as fo:
            hdul = fits.HDUList(hdus=[dust_densities_HDU] + additional_hdus)
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
        n_phi_new, n_rad_new = self.output_grid["rg"].shape
        assert n_rad_new, n_phi_new == self.io.OUT.shape[0, 2]

        density_keys = sorted(filter(lambda k: "rho" in k, self.input_data.fields.keys()))
        self._new_2D_arrays = np.array([self._interpolate2D(datakey=k) for k in density_keys])
        assert self._new_2D_arrays[0].shape == (n_rad_new, n_phi_new)

    @property
    def aspect_ratio(self):
        """Dimensionless ratio implied by mcfost parameters"""
        mcfl = self.config['mcfost_output']
        return mcfl['scale_height'] / mcfl['ref_radius']

    def gen_3D_arrays(self):
        '''Interpolate input data onto full 3D output grid'''
        nphi, nr = self.output_grid['rg'].shape
        nz_out, nr2 = self.output_grid['zg'].shape
        nz_in = self.config['mcfost_output']['nz']
        assert nr2 == nr
        assert nz_out == 2*nz_in + EXPECTED_ZSHAPE_INCREMENT

        nbins = len(self.new_2D_arrays)
        self._new_3D_arrays = np.zeros((nbins, nphi, nz_in, nr))
        for ir, r in enumerate(self.output_grid['rv']):
            z_vect = self.output_grid['zg'][nz_in+EXPECTED_ZSHAPE_INCREMENT:, ir].reshape(1, nz_in)
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

    cargs = parser.parse_args()

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
