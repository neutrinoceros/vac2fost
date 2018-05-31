'''Vizualize the app's output with pytest.'''

import os
import pathlib
import subprocess
import shutil

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import f90nml

from vac2fost import main as app

here = pathlib.Path(__file__).absolute().parent

# int tests -------------------------------------------------------------------

import pytest
@pytest.mark.incremental #each test is run only if the previous one passed
class TestLocalVersion():
    def test_python_call(self):
        app(str(here/'sample/vac2fost_conf.nml'))

    def test_format(self):
        f = pyfits.open(here.parent / 'hd142527_dusty0000.fits')[0]
        opt = f90nml.read(here / 'sample/vac2fost_conf.nml')['mcfost_list']
        assert f.data.shape[1:] == (opt['nphi'], 2*opt['nz']+1, opt['nr'])

@pytest.mark.skip(reason='unstable at the moment')
class TestInstalledVersion():
    '''require that vac2fost.py be accessible through your $PATH env variable'''
    def test_command_line_call(self):
        exitcode = subprocess.call(
            'vac2fost.py ' + str(here / 'sample/vac2fost_conf.nml'),
            shell=True
        )
        assert exitcode == 0

    def test_command_line_call_w_offset(self):
        exitcode = subprocess.call(
            'vac2fost.py ' + str(here / 'sample/vac2fost_conf.nml -o 2'),
            shell=True
        )
        assert exitcode == 0

# -----------------------------------------------------------------------------

if __name__=='__main__':
    '''Visualize the .fits data.

    Grid info is currently transmitted through app() outputs.
    '''
    out = app(str(here/'sample/vac2fost_conf.nml'))

    # get the Primary (only image available),
    # and exctract its first 3d array (density field)
    data = pyfits.open(out['finame'])[0].data[0]
    X = out['rads']*np.cos(out['phis'])
    Y = out['rads']*np.sin(out['phis'])
    fig, axes = plt.subplots(nrows=3, figsize=(10,15))

    nr, nz, nphi = data.shape
    vertical_profile = data[0,:,0]
    vertical_slice   = data[0,:,:]
    plane_slice      = data[:,int(nz/2),:]

    axes[0].plot(vertical_profile) #this should be a gaussian

    axes[1].set_aspect('equal')
    axes[1].pcolormesh(X, Y, plane_slice, cmap='inferno')

    axes[2].imshow(vertical_slice, cmap='inferno')

    fig.savefig('graph_test.png')

