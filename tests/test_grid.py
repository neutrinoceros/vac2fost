import subprocess
import pathlib
import shutil
import os

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from vac2fost import main as app


here = pathlib.Path(__file__).absolute().parent
output = here / 'output/'

import pytest
@pytest.mark.skip(reason='current version of mcfost does not detect absolute paths')
def test_grid_eq():
    res = app(str(here.parent / 'sample/vac2fost_conf.nml'))
    if pathlib.Path(output / 'data_disk').exists():
        shutil.rmtree(output / 'data_disk')
    subprocess.call(f'mcfost mcfost_conf.para -disk_struct -root_dir {output}', shell=True) #generate a grid data file
    grid = pyfits.open(str(output / 'data_disk/grid.fits.gz'))[0].data
    diffs = (res['rads'] - grid[0,:,0,:])[0]
    assert diffs.all() < 1e-13
