#!/usr/bin/env python3
'''A script for converting vac .dat data files to .fits

Steps
1) convert .dat to .blk with amrvac itself
2) convert .blk to .fits.gz
    a) sort the raw data points (original order is messy)
    b) reshape to 2D numpy arrays
    c) interpolate data to fit a MCFOST grid (log-spacing)
    d) convert to 3D (gaussian redistribution of density)
    e) unit conversion                                    <<<<< NOT IMPLEMENTED YET


Arguments
-c configuration file for this script
-d .dat file to be converted


Known limitations
   1) amr is not supported
   2) portability is not guaranted
   3) gaussian redistribution is currently unreliable (hard coded parameters)
   4) interpolation does not account for the curvature of polar cells
   5) a cylindrical grid is currently used for 3D,
      we may later implement the spherical option
   6) input simulation is assumed to be 2D (r,phi)
   7) gas density not being read yet
   8) when dust density is available, gas density is being ignored.
      That needs fixing if we wish to generate molecular lines synthetic observations.
'''

import pathlib
import subprocess
from argparse import ArgumentParser

import numpy as np
import astropy.io.fits as pyfits
from scipy.interpolate import interp2d
import f90nml

from amrvac_pywrap import interpret_shell_path, read_amrvac_conf
from vtk_vacreader import VacDataSorter
from vac2fost.utils import get_mcfost_grid

try:
    res = subprocess.check_output('which mcfost', shell=True).decode('utf-8')
    assert not 'not found' in res
except AssertionError:
    raise EnvironmentError('Installation of MCFOST not found.')


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

def get_grain_micron_sizes(amrvac_conf:f90nml.Namelist) -> np.ndarray:
    '''Read grain sizes (assumed in [cm]), from AMRVAC parameters and
    convert to microns.'''
    cm_sizes = np.array(amrvac_conf['usr_dust_list']['grain_size'])
    µm_sizes = 1e4 * cm_sizes
    return µm_sizes


def main(config_file:str, offset:int=None, output_dir:str='.', dbg=False):
    print('==========================================')
    print(          'Starting vac2fost.main()')
    print('==========================================')

    # .. input reading ..

    if isinstance(config_file, f90nml.Namelist):
        config = config_file
    else:
        config = f90nml.read(config_file)

    if offset is None:
        offset = config['target_options']['offset']
    outnum = str(offset).zfill(4)

    if isinstance(output_dir, str):
        output_dir = pathlib.Path(output_dir)

    options = config['target_options']
    sim_conf = read_amrvac_conf(files=options['amrvac_conf'], origin=options['origin'])
    vtu_filename = sim_conf['filelist']['base_filename'] + f'{outnum}.vtu'
    datfile = interpret_shell_path(options['origin']) + '/' + vtu_filename
    datshape = tuple([sim_conf['meshlist'][f'domain_nx{n}'] for n in (1,2)])

    print(f'loading data from {datfile}')
    simdata = VacDataSorter(file_name=datfile, data_shape=datshape)


    # .. interpolation to target grid ..

    target_grid = get_mcfost_grid(
        mcfost_list=config['mcfost_list'],
        mesh_list=sim_conf['meshlist'],
        silent=(not dbg)
    )
    rad_grid_new = target_grid[0,:,0,:].T
    phi_grid_new = target_grid[2,:,0,:].T
    n_rad_new, n_phi_new = rad_grid_new.shape
    assert n_rad_new == config['mcfost_list']['nr']
    assert n_phi_new == config['mcfost_list']['nphi']
    rad_vect_new = rad_grid_new[:,0]
    phi_vect_new = phi_grid_new[0]

    rad_vect_old, azim_vect_old = [simdata.get_axis(n) for n in range(2)]

    density_keys = sorted(filter(lambda k: 'rho' in k, simdata.fields.keys()))
    interpolated_arrays = []
    for k in density_keys:
        interpolator = interp2d(azim_vect_old, rad_vect_old, simdata[k], kind='cubic')
        interpolated_arrays.append(interpolator(phi_vect_new, rad_vect_new))
    assert interpolated_arrays[0].shape == (n_rad_new, n_phi_new)
    print("interpolation ok")


    # .. conversion to 3D ..

    zmax = config['target_options']['zmax']
    nz = config['mcfost_list']['nz']
    z_vect = np.linspace(-zmax, zmax, 2*nz+1)
    scale_height_grid = config['target_options']['aspect_ratio'] * rad_grid_new
    threeD_arrays = np.array([twoD2threeD(arr, scale_height_grid, z_vect) for arr in interpolated_arrays])
    print("3D conversion ok")


    # .. build a .fits file ..
    grain_sizes = get_grain_micron_sizes(sim_conf)
    assert len(grain_sizes) == len(threeD_arrays) - 1

    #the transposition is handling a weird behavior of fits files...
    dust_densities_array = np.stack(threeD_arrays[1 + grain_sizes.argsort()], axis=3).transpose()
    dust_densities_HDU = pyfits.PrimaryHDU(dust_densities_array)

    mcfost_keywords = {
        'read_n_a': 0, #automatic normalization of size-bins from mcfost param file.
        # following keywords are too long according to fits standards  !
        # --------------------------------------------------------------
        #'read_gas_density': 0, #set to 1 to add gas density
        #'gas_to_dust': sim.conf['usr_dust_list']['gas2dust_ratio'], #required when reading gas
    }

    for it in mcfost_keywords.items():
        dust_densities_HDU.header.append(it)

    grain_sizes_HDU = pyfits.ImageHDU(grain_sizes[grain_sizes.argsort()])

    hdus = [
        dust_densities_HDU,
        grain_sizes_HDU,
        #pyfits.ImageHDU(gas_density)
    ]
    fits_filename = output_dir / pathlib.Path(vtu_filename).name.replace('.vtu', '.fits')
    with open(fits_filename, 'w') as fo:
        hdul = pyfits.HDUList(hdus=hdus)
        hdul.writeto(fo)
    print(f'Successfully wrote {fits_filename}')

    print('==========================================')
    print(          'End of vac2fost.main()')
    print('==========================================')
    # .. finally, yield some info back (for testing) ..

    return dict(
        finame = fits_filename,
        rads   = rad_grid_new.T,
        phis   = phi_grid_new.T,
    )




if __name__=='__main__':
    # Parse the script arguments
    parser = ArgumentParser()
    parser.add_argument(
        dest='configuration', type=str,
        help='configuration file (namelist) for this script'
    )
    parser.add_argument(
        '-o', dest='off', type=str,
        required=False,
        default=None,
        help='output number of the target .dat VAC output file to be converted'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='activate debug mode (verbose)'
    )
    script_args = parser.parse_args()

    main(
        config_file=script_args.configuration,
        offset=script_args.off,
        dbg=script_args.debug
    )
