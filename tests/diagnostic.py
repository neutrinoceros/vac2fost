'''Visualize the .fits data.

Grid info is currently transmitted through app() outputs.
'''

import pathlib

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits

from vac2fost import main as app

here = pathlib.Path(__file__).absolute().parent

if __name__=='__main__':
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
