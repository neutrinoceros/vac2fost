'''Visualize the .fits data.

Grid info is currently transmitted through app() outputs.
'''

import pathlib

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits

from vac2fost import main as app

doc_dir = pathlib.Path(__file__).absolute().parent
conf = doc_dir.parent/"tests/sample/vac2fost_conf_nonaxisym.nml"
output_dir = doc_dir.parent/"demo_out"

if __name__=='__main__':
    itf = app(conf, output_dir=output_dir, verbose=False)

    # get the Primary (only image available),
    # and exctract its first 3d array (density field)
    filepath = itf.io['out'].filepath
    data = pyfits.open(filepath)[0].data[0]
    X = (itf.output_grid['rg'] * np.cos(itf.output_grid['phig'])).T
    Y = (itf.output_grid['rg'] * np.sin(itf.output_grid['phig'])).T

    nr, nz, nphi = data.shape
    vertical_profile = data[0,:,0]
    vertical_slice   = data[0,:,:]
    plane_slice      = data[:,int(nz/2),:]

    fig, axes = plt.subplots(nrows=3, figsize=(10,15))
    axes[0].set_title("mid plane slice (cartesian)")
    axes[0].set_aspect('equal')
    axes[0].pcolormesh(X, Y, plane_slice, cmap='inferno')

    axes[1].set_title("vertical distribution (gaussian)")
    axes[1].plot(vertical_profile, lw=0, marker="o")

    axes[2].set_title("vertical slice (not to scale)")
    axes[2].imshow(vertical_slice, cmap="inferno", origin="lower")

    fig.savefig(output_dir/'diagnostic_output.png')
