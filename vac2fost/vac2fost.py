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
  5) gas density is never passed to MCFOST as is but only
     as a tracer for smallest dust grains
'''
__version__ = '2.1.1'

from collections import OrderedDict as od, namedtuple
import os
from socket import gethostname
import sys
import shutil
import subprocess
from argparse import ArgumentParser
from pathlib import Path
import uuid

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp2d
import f90nml
try:
    import colorama
except ImportError:
    colorama = None

from amrvac_pywrap import interpret_shell_path, read_amrvac_conf
from vtk_vacreader import VacDataSorter

try:
    res = subprocess.check_output('which mcfost', shell=True).decode('utf-8')
    assert 'not found' not in res
except AssertionError:
    raise EnvironmentError('Installation of MCFOST not found.')

# Globals
MINGRAINSIZE_µ = 0.1
DEFAULTS = {'DBM': 'auto'}
DataInfo = namedtuple(
    'DataInfo',
    ['shape', 'directory', 'filename', 'filepath']
)

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
            ('Wavelengts', (
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
                od([('mol_data_file', 'co@xplot.dat'),
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

    def write_mcfost_conf(output_file: str, custom: dict = None, silent=True):
        '''Write a configuration file for mcfost using values from <custom>,
        and falling back to defaults found in block_descriptor defined above
        '''
        if custom is None:
            custom = {}
        if Path(output_file).exists() and not silent:
            print(f'Warning: {output_file} already exists, and will be overwritten.')
        with open(output_file, 'wt') as fi:
            fi.write('3.0'.ljust(10) + 'mcfost minimal version\n\n')
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
        if not silent:
            print(f'wrote {output_file}')

    def translate_amrvac_conf(itf) -> dict:
        # itf must be of type Interface (can't be parsed properly before python 3.7)
        '''pass amrvac parameters to mcfost'''
        parameters = {}

        # Zone
        mesh = itf.sim_conf['meshlist']
        parameters.update({
            'rin': mesh['xprobmin1']*itf.conv2au,
            'rout': mesh['xprobmax1']*itf.conv2au,
            'maps_size': 2*mesh['xprobmax1']*itf.conv2au,
        })

        try:
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
        except KeyError:
            itf.warnings.append("Could not find 'usr_dust_list' parameter, using default values")

        return parameters

    def get_mcfost_grid(itf) -> np.ndarray:
        '''Pre-run MCFOST with -disk_struct flag to get the exact grid used.'''
        mcfost_conf_file = itf.mcfost_para_file
        output_dir = itf.io['out'].directory
        silent = (not itf.dbg)

        output_dir = Path(output_dir).resolve()
        mcfost_conf_path = Path(mcfost_conf_file)
        if not output_dir.exists():
            subprocess.call(f'mkdir -p {output_dir}', shell=True)

        grid_file_name = output_dir / 'mcfost_grid.fits.gz'

        #gen_needed = True
        #mcfost_list = itf.config['mcfost_list']
        if grid_file_name.exists():
            itf.warnings.append("found existing grid file, ignored it")
            # devnote : this block is deprecated because it was getting off hand.
            #           a more maintainable solution is being studied.
            #
            # with fits.open(grid_file_name, mode='readonly') as fi:
            #     target_grid = fi[0].data
            # shape_found = target_grid.shape[1:]
            # correct_shapes = (
            #     (mcfost_list['nphi'], mcfost_list['nz'], mcfost_list['nr']),
            #     (mcfost_list['nphi'], mcfost_list['nz']*2+1, mcfost_list['nr'])
            # )
            # radial_range_correct = itf.conv2au*np.array([
            #     itf.sim_conf["meshlist"]["xprobmin1"],
            #     itf.sim_conf["meshlist"]["xprobmax1"]
            # ])
            # radial_range_found = np.array([target_grid[0, 0, 0, :].min(),
            #                                target_grid[0, 0, 0, :].max()])
            # zmax_found = target_grid[1, 0, :, :].max()
            # zmax_correct = itf.config["target_options"]["zmax"] * itf.conv2au
            # gen_needed = shape_found not in correct_shapes \
            #              or not np.all(radial_range_found == radial_range_correct) \
            #              or zmax_found != zmax_correct

        #if gen_needed: #devnote : always consider this true from now on
        assert mcfost_conf_path.exists()
        # generate a grid data file with mcfost itself and extract it
        tmp_mcfost_dir = Path(f'TMP_VAC2FOST_MCFOST_GRID_{uuid.uuid4()}')
        os.mkdir(tmp_mcfost_dir)
        try:
            shutil.copyfile(mcfost_conf_path.resolve(),
                            tmp_mcfost_dir/mcfost_conf_path.name)
        except shutil.SameFileError:
            pass

        pile = Path.cwd()
        os.chdir(tmp_mcfost_dir)
        try:
            os.environ['OMP_NUM_THREADS'] = '1'
            subprocess.check_call(
                f"mcfost mcfost_conf.para -disk_struct",
                shell=True,
                stdout={True: subprocess.PIPE, False: None}[silent]
            )
            shutil.move("data_disk/grid.fits.gz", grid_file_name)
        except subprocess.CalledProcessError as exc:
            errtip = f'\nError in MCFOST, exited with exitcode {exc.returncode}'
            if exc.returncode == 174:
                errtip += (
                    '\nThis is probably a memory issue. '
                    'Try reducing the target resolution or,'
                    ' alternatively, give more cpu memory to this task.'
                )
                raise RuntimeError(errtip)
        finally:
            os.chdir(pile)
            shutil.rmtree(tmp_mcfost_dir)
        with fits.open(grid_file_name, mode='readonly') as fi:
            target_grid = fi[0].data
        # end if gen_needed
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


def generate_conf_template() -> f90nml.Namelist:
    '''Generate a template namelist object with comments instead of default values'''
    target = {
        'origin': '!path to the simulation repository, where datafiles are located',
        'amrvac_conf': '!one or multiple file path relative to origin, ","separeated',
        'conv2au': 100,
        'num': 0
    }
    mcfost_list = {
        'nr': 128,
        'nr_in': 4,
        'nphi': 128,
        'nz': 10,
        # aspect ratio is implied by those parameters
        "flaring_index": 1.125,
        "ref_radius": 100.0,  # [a.u.]
        "scale_height": 1.0,  # [a.u.], at defined at ref_radius
    }
    template = f90nml.Namelist({
        'mcfost_list': f90nml.Namelist(mcfost_list),
        'target_options': f90nml.Namelist(target)
    })
    return template

# decorators
def parameterized(dec):
    """meta decorator, allow definition of decorators with parameters
    source: https://stackoverflow.com/questions/5929107/decorators-with-parameters"""
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)
        return repl
    return layer

@parameterized
def wait_for_ok(func, mess, lenght=61):
    """decorator, sandwich the function execution with '<mess>  ...' & 'ok'"""
    def modfunc(*args, **kwargs):
        print(mess.ljust(lenght), end="... ", flush=True)
        result = func(*args, **kwargs)
        print("ok")
        return result
    return modfunc

class Interface:
    '''A class to hold global variables as attributes and give
    clear and concise structure to the main() function.'''

    known_dbms = {'dust-only', 'gas-only', 'mixed', 'auto'}

    @wait_for_ok("parsing input")
    def __init__(self, config_file, num: int = None,
                 output_dir: Path = Path('.'),
                 dust_bin_mode: str = DEFAULTS['DBM'], dbg=False):

        # input checking
        if not isinstance(config_file, (str, Path)):
            raise TypeError(config_file)
        if not isinstance(output_dir, (str, Path)):
            raise TypeError(output_dir)
        if dust_bin_mode not in __class__.known_dbms:
            raise KeyError(f'Unknown dust binning mode "{dust_bin_mode}"')
        else:
            self._dust_binning_mode = dust_bin_mode

        # attribute storage
        self._base_args = {
            'config_file': Path(config_file),
            'output_dir': Path(output_dir),
            'num': num,
            'dust_bin_mode': dust_bin_mode,
        }

        self._dim = 2  # no support for 3D input yet
        self.dbg = dbg
        self.warnings = []

        # parse configuration file
        self.config = f90nml.read(config_file)
        if num is not None:
            self.num = num
        else:
            self.num = self.config['target_options']['num']

        origin = Path(self.config['target_options']['origin'])
        if not origin.is_absolute():
            to = self.config['target_options']
            p1 = Path.cwd()
            p2 = (Path(config_file).parent/origin).resolve()

            if isinstance(to['amrvac_conf'], (list, tuple)):
                fi = to['amrvac_conf'][0]
            else:
                fi = to['amrvac_conf']

            found = [(p/fi).is_file() for p in (p1, p2)]
            if all(found) and p1 != p2:
                raise FileNotFoundError(
                    f"""can not guess if path "{origin}" is relative to cwd or configuration file""")
            elif not any(found):
                raise FileNotFoundError(origin)
            else:
                p = (p1, p2)[found.index(True)]
            self.config['target_options'].update({'origin': p.resolve()})
        self.sim_conf = read_amrvac_conf(
            files=self.config['target_options']['amrvac_conf'],
            origin=self.config['target_options']['origin']
        )
        self.small_grains_from_gas = True
        self._iodat = None
        self._input_data = None
        self._output_grid = None
        self._µsizes = None
        self._new_2D_arrays = None
        self._new_3D_arrays = None

        if not self.io['out'].directory.exists():
            subprocess.call(f"mkdir -p {self.io['out'].directory}",
                            shell=True)
            self.warnings.append(f"rep {self.io['out'].directory} was created")

        # optional definition of the distance unit
        self.conv2au = 1.0
        try:
            self.conv2au = self.config['target_options']['conv2au']
        except KeyError:
            self.warnings.append("could not find conv2au, distance unit assumed 1au")

    def print_warnings(self):
        '''Print warnings if any.'''
        if self.warnings:
            if colorama is not None:
                colorama.init()
                red = colorama.Fore.RED + colorama.Style.BRIGHT
            else: red = ""
            print()
            print(" WARNINGS:")
            print(red+'\n'.join([f" - {w}" for w in self.warnings]))
            if colorama is not None:
                print(colorama.Style.RESET_ALL, end='')

    @property
    def dust_binning_mode(self):
        return self._dust_binning_mode

    @dust_binning_mode.setter
    def dust_binning_mode(self, args:list):
        if isinstance(args, str):
            dbm = args
            reason = None
        else:
            assert isinstance(args, list) and len(args) == 2
            dbm, reason = args
        assert dbm in __class__.known_dbms
        w = f'dust-binning mode was switched to "{args[0]}"'
        if reason:
            w += f'; reason: {reason}'
        self.warnings.append(w)
        self._dust_binning_mode = dbm

    @property
    def grain_micron_sizes(self) -> np.ndarray:
        '''Read grain sizes (assumed in [cm]), from AMRVAC parameters and
        convert to microns.'''
        µm_sizes = np.empty(0)
        if self._µsizes is None:
            if self.dust_binning_mode != 'gas-only':
                try:
                    cm_sizes = np.array(
                        self.sim_conf['usr_dust_list']['grain_size_cm'])
                    µm_sizes = 1e4 * cm_sizes
                except KeyError:
                    if self.dust_binning_mode == "auto":
                        self.dust_binning_mode = ["gas-only", "could not find grain sizes"]
                    else:
                        raise

            if min(µm_sizes) > MINGRAINSIZE_µ and self.dust_binning_mode == "auto":
                self.dust_binning_mode = ["mixed",
                                          f"smallest size found > {MINGRAINSIZE_µ}µm"]

            if self.dust_binning_mode in {'gas-only', 'mixed'}:
                µm_sizes = np.insert(µm_sizes, 0, MINGRAINSIZE_µ)
            self._µsizes = µm_sizes
        return self._µsizes

    @property
    def argsort_offset(self):
        '''Get the slice starting index when selecting arrays to be transformed'''
        return 1 - int(self.small_grains_from_gas)

    @property
    def io(self) -> dict:
        '''Store general info on input/output file locations
        and data array shapes.'''
        if self._iodat is None:
            vtu_filename = ''.join([self.sim_conf['filelist']['base_filename'],
                                    str(self.num).zfill(4),
                                    '.vtu'])
            self._iodat = {}
            basein = dict(
                directory=Path(interpret_shell_path(self.config['target_options']['origin'])).resolve(),
                filename=vtu_filename,
                shape=tuple(
                    [self.sim_conf['meshlist'][f'domain_nx{n}']
                     for n in range(1, self._dim+1)]
                )
            )
            baseout = dict(
                directory=Path(self._base_args['output_dir']),
                filename=basein['filename'].replace('.vtu', '.fits'),
                shape=None  # not used: don't write bugs when you don't need to
            )
            for d, k in zip([basein, baseout], ['in', 'out']):
                d.update(
                    {'filepath': (d['directory'] / d['filename']).resolve()})
                self._iodat.update({k: DataInfo(**d)})
        return self._iodat

    @property
    def mcfost_para_file(self):
        '''Locate output configuration file for mcfost'''
        file = self.io['out'].directory/'mcfost_conf.para'
        return str(file)

    def load_input_data(self) -> None:
        '''Use vtkvacreader.VacDataSorter to load AMRVAC data'''
        self._input_data = VacDataSorter(
            file_name=str(self.io['in'].filepath),
            shape=self.io['in'].shape
        )

    @property
    def input_data(self):
        '''Load input simulation data'''
        if self._input_data is None:
            self.load_input_data()
        return self._input_data

    @property
    def output_grid(self) -> dict:
        '''Store info on 3D output grid specifications
        as vectors "v", and (r-phi)grids "g"'''
        if self._output_grid is None:
            if not Path(self.mcfost_para_file).is_file():
                self.write_mcfost_conf_file()
            target_grid = MCFOSTUtils.get_mcfost_grid(self)
            self._output_grid = {
                'array': target_grid,
                # (nr, nphi) 2D grids
                'rg': target_grid[0, :, 0, :].T,
                'phig': target_grid[2, :, 0, :].T,
                # (nr, nz) 2D grid (z points do not depend on phi)
                'zg': target_grid[1, 0, :, :].T,
                # vectors (1D arrays)
                'rv': target_grid[0, :, 0, :].T[:, 0],
                'phiv': target_grid[2, :, 0, :].T[0],
            }
        return self._output_grid

    def write_mcfost_conf_file(self) -> None:
        '''Customize defaults with user specifications'''
        custom = {}
        custom.update(MCFOSTUtils.translate_amrvac_conf(self))
        unknown_args = self.scan_for_unknown_arguments()
        if unknown_args:
            raise KeyError(f'Unrecognized MCFOST argument(s): {unknown_args}')
        custom.update(self.config['mcfost_list'])

        custom.update({'dust_mass': get_dust_mass(self.input_data)})
        MCFOSTUtils.write_mcfost_conf(
            output_file=self.mcfost_para_file,
            custom=custom,
            silent=(not self.dbg)
        )

    def scan_for_unknown_arguments(self) -> list:
        '''Get unrecognized arguments found in mcfost_list'''
        unknowns = []
        for arg in self.config['mcfost_list'].keys():
            if not arg.lower() in MCFOSTUtils.known_args:
                unknowns.append(arg)
        return unknowns

    def write_output(self) -> None:
        '''Main method. Write a .fits file suited for MCFOST input.'''
        # the transposition is handling a weird behavior of fits files...
        dust_densities_array = np.stack(
            self.new_3D_arrays[
                self.argsort_offset + self.grain_micron_sizes.argsort()],
            axis=3).transpose()
        dust_densities_HDU = fits.PrimaryHDU(dust_densities_array)

        mcfost_keywords = {
            # automatic normalization of size-bins from mcfost param file.
            'read_n_a': 0,
            # following keywords are too long according to fits standards !
            # issue 19
            # -------------------------------------------------------------
            # 'read_gas_density': 0, #set to 1 to add gas density
            # required when reading gas
            # 'gas_to_dust': sim.conf['usr_dust_list']['gas2dust_ratio'],
        }

        for it in mcfost_keywords.items():
            dust_densities_HDU.header.append(it)

        grain_sizes_HDU = fits.ImageHDU(
            self.grain_micron_sizes[self.grain_micron_sizes.argsort()]
        )
        hdus = [
            dust_densities_HDU,
            grain_sizes_HDU,
            # fits.ImageHDU(gas_density) # issue 19 related...
        ]
        fopath = self.io['out'].filepath
        with open(fopath, 'wb') as fo:
            hdul = fits.HDUList(hdus=hdus)
            hdul.writeto(fo)

    @property
    def input_grid(self) -> dict:
        '''Store physical coordinates (vectors)
        about the input grid specifications.'''
        ig = {
            'rv': self.input_data.get_ticks('r') * self.conv2au,
            'phiv': self.input_data.get_ticks('phi')
        }
        return ig

    def gen_2D_arrays(self):
        '''Interpolate input data onto r-phi grid
        with output grid specifications'''
        n_rad_new, n_phi_new = self.output_grid['rg'].shape
        assert n_rad_new == self.config['mcfost_list']['nr']
        assert n_phi_new == self.config['mcfost_list']['nphi']

        density_keys = sorted(filter(
            lambda k: 'rho' in k, self.input_data.fields.keys()))
        interpolated_arrays = []
        for k in density_keys:
            interpolator = interp2d(
                self.input_grid['phiv'],
                self.input_grid['rv'],
                self.input_data[k], kind='cubic'
            )
            interpolated_arrays.append(
                interpolator(self.output_grid['phiv'],
                             self.output_grid['rv'])
            )
        assert interpolated_arrays[0].shape == (n_rad_new, n_phi_new)
        self._new_2D_arrays = np.array(interpolated_arrays)

    @property
    def aspect_ratio(self):
        """Dimensionless ratio implied by mcfost parameters"""
        mcfl = self.config['mcfost_list']
        return mcfl['scale_height'] / mcfl['ref_radius']

    def gen_3D_arrays(self):
        '''Interpolate input data onto full 3D output grid'''
        nr, nphi = self.output_grid['rg'].shape
        nr2, nz_out = self.output_grid['zg'].shape
        nz_in = self.config['mcfost_list']['nz']
        assert nr2 == nr
        assert nz_out == 2*nz_in+1

        nbins = len(self.new_2D_arrays)
        self._new_3D_arrays = np.zeros((nbins, nr, nz_in, nphi))
        for ir, r in enumerate(self.output_grid['rv']):
            z_vect = self.output_grid['zg'][ir, nz_in+1:].reshape(1, nz_in)
            local_height = r * self.aspect_ratio
            gaussian = np.exp(-z_vect**2/ (2*local_height**2)) / (np.sqrt(2*np.pi) * local_height)
            for i_bin, surface_density in enumerate(self.new_2D_arrays[:, ir, :]):
                self._new_3D_arrays[i_bin, ir, :, :] = np.transpose(
                    gaussian * surface_density.reshape(nphi, 1))

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

# =======================================================================================
class VerbatimInterface(Interface):
    """A more talkative Interface"""
    @wait_for_ok(f"loading input data")
    def load_input_data(self) -> None:
        super().load_input_data()

    @wait_for_ok('writting mcfost configuration file')
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


# =======================================================================================
def main(config_file: str,
         num: int = None,
         output_dir: str = '.',
         dust_bin_mode: str = DEFAULTS['DBM'],
         verbose=False,
         dbg=False):
    '''Try to transform a .vtu file into a .fits'''

    print('=========================== vac2fost.py ============================')
    InterfaceType = {True: VerbatimInterface, False: Interface}[verbose]
    itf = InterfaceType(config_file, num=num, output_dir=output_dir,
                        dust_bin_mode=dust_bin_mode, dbg=dbg)

    itf.load_input_data()
    itf.write_mcfost_conf_file()
    itf.gen_2D_arrays()
    itf.gen_3D_arrays()
    itf.write_output()
    itf.print_warnings()
    print(f"\nsuccess ! output wrote:\n{itf.io['out'].filepath}")

    print('=========================== end program ============================')

    # return the Interface object for inspection (tests)
    return itf


# =======================================================================================
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
        '-n', dest='num', type=int,
        required=False,
        default=None,
        help='output number of the target .vtu VAC output file to be converted'
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
        default=DEFAULTS['DBM'],
        help='prefered bin selection mode (accepted values "dust-only", "gas-only", "mixed", "auto")'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='activate verbose mode'
    )
    parser.add_argument(
        '--dbg', '--debug', dest='dbg',
        action='store_true',
        help='activate debug mode (verbose for MCFOST)'
    )
    parser.add_argument(
        '--genconf', action='store_true',
        help="print a default configuration file for vac2fost"
    )
    parser.add_argument(
        '--profile',
        action='store_true',
        help='activate profiling mode'
    )

    cargs = parser.parse_args()

    if cargs.genconf:
        print(generate_conf_template())
        print(f"%% automatically generated with vac2fost {__version__}\n")
        sys.exit(0)
    elif len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if cargs.profile:
        import cProfile
        import pstats
        import io
        pr = cProfile.Profile()
        pr.enable()
    # -------------------------------------------
    main(
        config_file=cargs.configuration,
        num=cargs.num,
        output_dir=cargs.output,
        dust_bin_mode=cargs.dbm,
        verbose=cargs.verbose,
        dbg=cargs.dbg
    )
    # -------------------------------------------
    if cargs.profile:
        pr.disable()
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
        ps.print_stats()
        print(s.getvalue())
