#!/usr/bin/env python3
'''A conversion tool from AMRVAC (.vtu) to MCFOST (.fits)

Run `vac2fost.py --help` for documentation on command line usage

Steps
    a) load AMRVAC data with vtk_vacreader.VacDataSorter(), sort it as 2D arrays
    b) dry-run MCFOST to get the exact target grid
    c) interpolate data to the target grid (log-spacing)
    d) convert to 3D (gaussian redistribution of density)

Known limitations
   1) amr is not supported by reader
   2) portability is not guaranted
   3) interpolation does not account for the curvature of polar cells
   4) a cylindrical grid is currently used for 3D,
      we may later implement the spherical option
   5) input simulation is assumed to be 2D (r,phi)
   6) gas density not being read yet
   7) when dust density is available, gas density is being ignored.
      That needs fixing if we wish to generate molecular lines synthetic observations.
'''

from collections import OrderedDict as od, namedtuple
import os
import sys
import shutil
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp2d
import f90nml

from amrvac_pywrap import interpret_shell_path, read_amrvac_conf
from vtk_vacreader import VacDataSorter

try:
    res = subprocess.check_output('which mcfost', shell=True).decode('utf-8')
    assert not 'not found' in res
except AssertionError:
    raise EnvironmentError('Installation of MCFOST not found.')

minsize_grain_µm = 0.1

class MCFOSTUtils:
    '''Utility functions to call MCFOST in vac2fost.main() to define the final grid.'''

    blocks_descriptors = od(
        #name every mcfost parameter (by order of appearence) and give default values
        [
            ('Photons', (
                od([('nphot_temp', '1.28e5')]),
                od([('nphot_sed', '1.28e3')]),
                od([('nphot_img', '1.28e5')])
            )),

            ('Wavelengts', (
                od([('n_lambda', 50), ('lambda_min', 0.1), ('lambda_max', 3e3)]),
                od([('compute_temp', True), ('compute_sed', True), ('use_default_wl', True)]),
                od([('wavelength_file', 'wavelengths.dat')]),
                od([('separation', False), ('stokes_parameters', False)])
            )),
            ('Grid', (
                 od([('geometry', '1')]),
                 od([('nr', 100), ('nz', 10), ('nphi', 100), ('nr_in', 30)])
            )),
            ('Maps', (
                od([('nx', 501), ('ny', 501), ('maps_size', 400)]),
                od([('imin', 0), ('imax', 0), ('n_incl', 1), ('centered', False)]),
                od([('az_min', 0), ('az_max', 240), ('n_az_angles', 1)]),
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
                od([('dust_settling', 3), ('exp_strat', 0.5), ('a_srat', 1.0)]),
                od([('dust_radial_migration', False)]),
                od([('sublimate_dust', False)]),
                od([('hydrostatic_eq', False)]),
                od([('viscous_heating', False), ('alpha_viscosity', '1e-3')]),
             )),
            ('Number of zones', (
                od([('n_zones', '1')]),
            )),
            ('Zone', (
                od([('zone_type', 1)]),
                od([('dust_mass', '1e-3'), ('gas_to_dust_ratio', 100)]),
                od([('scale_height', 10.0), ('ref_radius', 100.0), ('profile_exp', 2)]),
                od([('rin', 10), ('edge', 0), ('rout', 200), ('rc', 100)]),
                od([('flaring_index', 1.125)]),
                od([('density_exp', -0.5), ('gamma_exp', 0.0)])
            )),
            ('Grains', (
                od([('n_species', 1)]),
                od([('grain_type', 'Mie'), ('n_components', 1), ('mixing_rule', 2), ('porosity', 0.), ('mass_fraction', 0.75), ('vmax_dhs', 0.9)]),
                od([('optical_indices_file', 'Draine_Si_sUV.dat'), ('volume_fraction', 1.0)]),
                od([('heating_method', 1)]),
                od([('sp_min', minsize_grain_µm), ('sp_max', 1000), ('sexp', 3.5), ('n_grains', 100)])
            )),
            ('Molecular RT', (
                od([('lpop', True), ('laccurate_pop', True), ('LTE', True), ('profile_width', 15.)]),
                od([('v_turb', 0.2)]),
                od([('nmol', 1)]),
                od([('mol_data_file', 'co@xplot.dat'), ('level_max', 6)]),
                od([('vmax', 1.0), ('n_speed', 20)]),
                od([('cst_mol_abund', True), ('abund', '1e-6'), ('abund_file', 'abundance.fits.gz')]),
                od([('ray_tracing', True), ('n_lines_rt', 3)]),
                od([('transition_num_1', 1), ('transition_num_2', 2), ('transition_num_3', 3)])
            )),
            ('Star', (
                od([('n_stars', 1)]),
                od([('star_temp', 4000.0), ('star_radius', 2.0), ('star_mass', 1.0), ('star_x',0.), ('star_y', 0.), ('star_z', 0), ('star_is_bb', True)]),
                od([('star_rad_file', 'lte4000-3.5.NextGen.fits.gz')]),
                od([('fUV', 0.0), ('slope_fUV', 2.2)]),
            ))
        ])

    def write_mcfost_conf(output_file:str, custom:dict={}, silent=True):
        '''Write a configuration file for mcfost using values from <custom>,
        and falling back to defaults found in block_descriptor defined above
        '''
        if Path(output_file).exists() and not silent:
            print(f'Warning: {output_file} already exists, and will be overwritten.')
        with open(output_file, 'w') as fi:
            fi.write('3.0'.ljust(10) + 'mcfost minimal version' + '\n\n')
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
                    fi.write('  ' + '  '.join(parameters).ljust(36) + '  ' + ', '.join(line.keys()))
                    fi.write('\n')
                fi.write('\n')
            fi.write(f'\n\n\n%% GENERATED BY {__file__} %%\n')
        if not silent:
            print(f'wrote {output_file}')

    def translate_amrvac_conf(itf) -> dict:
        #itf must be of type Interface (can't be parsed properly before python 3.7
        '''pass amrvac parameters to mcfost'''
        parameters = {}

        # Zone
        mesh = itf.sim_conf['meshlist']
        parameters.update({
            'rin': mesh['xprobmin1']*itf.conv2au,
            'rout': mesh['xprobmax1']*itf.conv2au,
            'maps_size': 2*mesh['xprobmax1']*itf.conv2au,
        })

        # aspect ratio may be defined in the hd simulation conf file
        try:
            parameters.update({
                'ref_radius': 1.0, #AU
                'scale_height': itf.sim_conf['disk_list']['aspect_ratio'] #at ref radius
            })
        except KeyError:
            pass

        try:
            dl2 = itf.sim_conf['usr_dust_list']
            parameters.update({
                'gas_to_dust_ratio': dl2['gas2dust_ratio'],
                #'dust_mass': ... #can not be passed from the configuration file alone
            })
            # Grains
            sizes_µm = itf.grain_micron_sizes
            parameters.update({
                #min/max grain sizes in microns
                'sp_min': min(1e-1, min(sizes_µm)),
                'sp_max': max(1e3,  max(sizes_µm)),
            })
        except KeyError:
            #in case the list 'usr_dust_list' is not found, pass default values to mcfost
            pass

        return parameters


    def get_mcfost_grid(mcfost_conf:str, mcfost_list:dict={}, output_dir:str='.', silent=True) -> np.ndarray:
        '''pre-run MCFOST in -disk_struct mode to extract the exact grid used.'''
        output_dir = Path(output_dir)
        if not output_dir.exists():
            subprocess.call(f'mkdir --parents {output_dir}', shell=True)

        grid_file_name = Path(output_dir) / 'mcfost_grid.fits.gz'

        gen_needed = True
        if grid_file_name.exists():
            with fits.open(grid_file_name, mode='readonly') as fi:
                target_grid = fi[0].data
            found = target_grid.shape
            hoped = mcfost_list['nphi'], mcfost_list['nz'], mcfost_list['nr']
            gen_needed = found[1:] != hoped

        if gen_needed:
            assert Path(mcfost_conf).exists()
            try:
                shutil.copyfile(output_dir / 'mcfost_conf.para', './mcfost_conf.para')
            except shutil.SameFileError:
                pass

            # generate a grid data file with mcfost itself and extract it
            tmp_fost_dir = Path('TMP_VAC2FOST_MCFOST_GRID')
            try:
                os.environ['OMP_NUM_THREADS'] = '1'
                subprocess.check_call(
                    f'mcfost mcfost_conf.para -disk_struct -root_dir {tmp_fost_dir}',
                    shell=True,
                    stdout={True: subprocess.PIPE, False: None}[silent]
                )
                if tmp_fost_dir.exists():
                    shutil.move(tmp_fost_dir / 'data_disk/grid.fits.gz', grid_file_name)
            except subprocess.CalledProcessError as exc:
                errtip = f'\nError in mcfost, exited with exitcode {exc.returncode}'
                if exc.returncode == 174:
                    errtip += (
                        '\nThis is probably a memory issue. '
                        'Try reducing your target resolution or alternatively, '
                        'give more cpu memory to this task.'
                    )
                raise RuntimeError(errtip)
            finally:
                if output_dir != Path('.'):
                    os.remove('./mcfost_conf.para')
                if tmp_fost_dir.exists():
                    shutil.rmtree(tmp_fost_dir)
            with fits.open(grid_file_name, mode='readonly') as fi:
                target_grid = fi[0].data
        return target_grid


def gauss(z, sigma):
    return 1./(np.sqrt(2*np.pi) * sigma) * np.exp(-z**2/(2*sigma**2))

def twoD2threeD(arr2d:np.ndarray, scale_height:np.ndarray, zvect:np.ndarray) -> np.ndarray:
    '''Convert surface density 2d array into volumic density 3d
    cylindrical array assuming a gaussian vertical distribution.

    formats
    arr2d : (nr, nphi)
    arr3d : (nr, nz, nphi) (suited for mcfost)

    note
    MCFOST offers the possibility to use a spherical grid instead.
    '''
    #devnote : gaussian distribution of dust is a bad fit.
    #For better modelization, see
    #eq 1 from (Pinte et al 2008) and eq 25 from (Fromang & Nelson 2009)

    nrad, nphi = arr2d.shape
    nz = len(zvect)
    arr3d = np.ones((nrad, nz, nphi))

    for k,z in enumerate(zvect):
        arr3d[:,k,:] = arr2d[:,:] * gauss(z, sigma=scale_height)
    return arr3d


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
        mass += np.sum([cell_surfaces * field[:,i] for i in range(field.shape[1])])
    return mass

def main(
        config_file:str,
        offset:int=None,
        output_dir:str='.',
        g2d_bin=False,
        read_gas=False,
        verbose=False,
        dbg=False
):
    printer = {
        True: print,
        False: lambda *args, **kwargs: None,
    }[verbose]


    printer(' --------- Start vac2fost.main() ---------')
    printer('reading input ...', end=' ', flush=True)
    itf = Interface(config_file, num=offset, output_dir=output_dir, g2d_bin=g2d_bin)
    printer('ok')

    printer(f"loading data from {itf.io['in'].filename}", end=' ', flush=True)
    simdata = itf.input_data
    printer('ok')

    printer('writting the mcfost configuration file ...', end=' ', flush=True)
    itf.write_mcfost_conf_file()
    printer('ok')

    # -------------------------------------------------------------

    printer('interpolating to MCFOST grid ...', end=' ', flush=True)
    target_grid = MCFOSTUtils.get_mcfost_grid(
        mcfost_conf=itf.mcfost_para_file,
        mcfost_list=itf.config['mcfost_list'],
        output_dir=itf.io['out'].directory,
        silent=(not dbg)
    )
    rad_grid_new = target_grid[0,:,0,:].T
    phi_grid_new = target_grid[2,:,0,:].T
    n_rad_new, n_phi_new = rad_grid_new.shape
    assert n_rad_new == itf.config['mcfost_list']['nr']
    assert n_phi_new == itf.config['mcfost_list']['nphi']
    rad_vect_new = rad_grid_new[:,0]
    phi_vect_new = phi_grid_new[0]

    rad_vect_old  = simdata.get_ticks('r') * itf.conv2au
    azim_vect_old = simdata.get_ticks('phi')

    density_keys = sorted(filter(lambda k: 'rho' in k, simdata.fields.keys())) #todo : update me
    interpolated_arrays = []
    for k in density_keys:
        interpolator = interp2d(azim_vect_old, rad_vect_old, simdata[k], kind='cubic')
        interpolated_arrays.append(interpolator(phi_vect_new, rad_vect_new))
    assert interpolated_arrays[0].shape == (n_rad_new, n_phi_new)
    printer('ok')


    # -------------------------------------------------------------
    printer('converting 2D arrays to 3D ...', end=' ', flush=True)
    zmax = itf.config['target_options']['zmax']
    nz = itf.config['mcfost_list']['nz']
    z_vect = np.linspace(0, zmax, nz)
    scale_height_grid = itf.config['target_options']['aspect_ratio'] * rad_grid_new
    threeD_arrays = np.array([twoD2threeD(arr, scale_height_grid, z_vect) for arr in interpolated_arrays])
    printer('ok')



    # -------------------------------------------------------------
    printer('building the .fits file ...', end=' ', flush=True)

    #the transposition is handling a weird behavior of fits files...
    dust_densities_array = np.stack(threeD_arrays[itf.argsort_offset + itf.grain_micron_sizes.argsort()], axis=3).transpose()
    dust_densities_HDU = fits.PrimaryHDU(dust_densities_array)

    mcfost_keywords = {
        'read_n_a': 0, #automatic normalization of size-bins from mcfost param file.
        # following keywords are too long according to fits standards  !
        # --------------------------------------------------------------
        #'read_gas_density': 0, #set to 1 to add gas density
        #'gas_to_dust': sim.conf['usr_dust_list']['gas2dust_ratio'], #required when reading gas
    }

    for it in mcfost_keywords.items():
        dust_densities_HDU.header.append(it)

    grain_sizes_HDU = fits.ImageHDU(itf.grain_micron_sizes[itf.grain_micron_sizes.argsort()])

    hdus = [
        dust_densities_HDU,
        grain_sizes_HDU,
        #fits.ImageHDU(gas_density)
    ]

    with open(itf.io['out'].filename, 'wb') as fo:
        hdul = fits.HDUList(hdus=hdus)
        hdul.writeto(fo)
    printer('ok')
    printer(f"Successfully wrote {itf.io['out'].filename}")
    printer(' --------- End   vac2fost.main() ---------')

    printer('Messages collection:')
    printer('\n    '.join(itf.messages))
    printer()
    printer('Warnings collection:')
    printer('\n    '.join(itf.warnings))

    # .. finally, yield some info back (for testing) ..
    return dict(
        finame = itf.io['out'].filename,
        rads   = rad_grid_new.T,
        phis   = phi_grid_new.T,
    )

def generate_conf_template():
    target = {
        'origin': '<path to the simulation repository, where datafiles are located>',
        'amrvac_conf': '<one or multiple file path relative to origin, separated by comas ",">',
        'zmax': '<<real> maximum height of the disk for vertical extrapolation (cylindrical). Use same unit as in .vtu data>',
        'aspect_ratio': '<<real> constant aspect ratio for vertical extrapolation>'
    }
    mcfost_params = {
        'nr': 128,
        'nr_in': 4,
        'nphi': 128,
        'nz': 10
    }
    template = f90nml.Namelist({
        'mcfost_list': f90nml.Namelist(mcfost_params),
        'target_options': f90nml.Namelist(target)
    })
    return template

#////////////////////////////////////////////////////////////
# v2.0 dev zone
DataInfo = namedtuple('DataInfo', ['shape', 'directory', 'filename'])

class Interface:
    '''a class to hold global variables as attributes and give
    better structure to the sequence'''

    def __init__(self, config_file:str, num:int=None, output_dir:str='.',
                 g2d_bin=False, dbg=False):
        self._base_args = {
            'config_file': config_file,
            'output_dir': output_dir,
            'num': num,
            'g2d_bin': g2d_bin,
            'dbg': dbg,
        }
        self._dim = 2 #no support for 3D input yet
        self.messages = []
        self.warnings = []

        if isinstance(config_file, f90nml.Namelist):
            self.config = config_file
        else:
            self.config = f90nml.read(config_file)

        self.num = num or self.config['target_options']['offset']

        to = self.config['target_options'] # alias
        self.sim_conf = read_amrvac_conf(files=to['amrvac_conf'], origin=to['origin'])

        self.small_grains_from_gas = True
        self._input_data = None

        if not self.io['out'].directory.exists():
            subprocess.call(f"mkdir --parents self.io['out'].directory", shell=True)
            self.warnings.append(f"rep {self.io['out'].directory} was created")

        if not (self.io['in'].directory/self.io['in'].filename).exists():
            raise FileNotFoundError(self.io['in'].directory/self.io['in'].filename)

        #optional definition of the distance unit
        self.conv2au = 1.0
        try:
            self.conv2au = self.config['target_options']['conv2au']
        except KeyError:
            self.warnings.append('parameter conv2au was not found. Distance unit in simulation is assumed to be 1au (astronomical unit).')

    @property
    def grain_micron_sizes(self) -> np.ndarray:
        '''Read grain sizes (assumed in [cm]), from AMRVAC parameters and
        convert to microns.'''
        try:
            cm_sizes = np.array(self.sim_conf['usr_dust_list']['grain_size_cm'])
            µm_sizes = 1e4 * cm_sizes
        except KeyError:
            #no grain size detected, gas only
            µm_sizes = []

        # decide if an additional fake dust bin, based on gas density, is necessary
        if µm_sizes != [] :
            if (min(µm_sizes) > minsize_grain_µm or self.g2d_bin):
                self.small_grains_from_gas = True
        else:
            self.small_grains_from_gas = True
            self.warnings.append('no grain size detected, using gas as a proxy')

        if self.small_grains_from_gas:
            µm_sizes = np.insert(µm_sizes, 0, minsize_grain_µm)
        #assert len(self.grain_micron_sizes) == len(self.threeD_arrays) - self.argsort_offset
        return µm_sizes

    @property
    def argsort_offset(self):
        return 1 - int(self.small_grains_from_gas)

    @property
    def io(self) -> dict:
        vtu_filename = self.sim_conf['filelist']['base_filename'] + str(self.num).zfill(4) + '.vtu'
        res = {}
        res.update({
            'in': DataInfo(
                directory=Path(interpret_shell_path(self.config['target_options']['origin'])),
                filename=vtu_filename,
                shape=tuple(
                    [self.sim_conf['meshlist'][f'domain_nx{n}'] for n in range(1, self._dim+1)]
                )
            )
        })
        res.update({
            'out': DataInfo(
                directory=Path(self._base_args['output_dir']),
                filename=res['in'].filename.replace('.vtu', '.fits'),
                shape='' #TODO : fill me
            )
        })
        return res

    @property
    def mcfost_para_file(self):
        return str(self.io['out'].directory/'mcfost_conf.para')

    @property
    def input_data(self):
        if self._input_data is None:
            self._input_data = VacDataSorter(
                file_name=str(self.io['in'].directory/self.io['in'].filename),
                shape=self.io['in'].shape
            )
        return self._input_data

    def write_mcfost_conf_file(self):
        custom = {}
        custom.update(MCFOSTUtils.translate_amrvac_conf(self))
        custom.update(self.config['mcfost_list'])
        custom.update({'dust_mass': get_dust_mass(self.input_data)})
        MCFOSTUtils.write_mcfost_conf(
            output_file=self.mcfost_para_file,
            custom=custom,
            silent=(not self._base_args['dbg'])
        )

    def run(self):
        pass
#////////////////////////////////////////////////////////////


if __name__=='__main__':
    # Parse the script arguments
    p = ArgumentParser(description='Parse arguments for main app')
    p.add_argument(
        dest='configuration', type=str,
        nargs='?',
        default=None,
        help='configuration file (namelist) for this script'
    )
    p.add_argument(
        '-n', dest='num', type=int,
        required=False,
        default=None,
        help='output number of the target .vtu VAC output file to be converted'
    )
    p.add_argument(
        '-o', '--output', dest='output', type=str,
        required=False,
        default='.',
        help='select output directory for generated files'
    )
    p.add_argument(
        '--g2d',
        action='store_true',
        help='activate gas-to-dust mode'
    )
    p.add_argument(
        '--gas',
        action='store_true',
        help='pass information on gas component to mcfost (not implemented !)'
    )
    p.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='activate verbose mode'
    )
    p.add_argument(
        '--dbg', '--debug', dest='dbg',
        action='store_true',
        help='activate debug mode (verbose for MCSOST)'
    )
    p.add_argument(
        '--genconf', action='store_true',
        help='generate configuration file template for this script in the current dir'
    )

    args = p.parse_args()

    if args.genconf:
        template = generate_conf_template()
        finame = args.output + '/template_vac2fost.nml'
        if not Path(args.output).exists():
            subprocess.call(f'mkdir --parents {args.output}', shell=True)
        if Path(finame).exists():
            sys.exit(f'Error: {finame} already exists, exiting vac2fost.py')
        else:
            with open(finame, 'w') as fi:
                template.write(fi)
                print(f'Generated {finame}')
        sys.exit()
    elif not args.configuration:
        sys.exit('Error: a configuration file is required as first argument. You can generate a template with --genconf')

    main(
        config_file=args.configuration,
        offset=args.num,
        output_dir=args.output,
        g2d_bin=args.g2d,
        read_gas=args.gas,
        verbose=args.verbose,
        dbg=args.dbg
    )
