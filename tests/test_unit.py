import pathlib
import pickle
import numpy as np
from astropy.io import fits as pyfits

from vac2fost.vac2fost import twoD2threeD
test_dir = pathlib.Path(__file__).resolve().parent

NR, NPHI, NZ = 3, 5, 7
rvect = np.linspace(1, 10, NR)
phivect = np.linspace(0, 2*np.pi, NPHI)
zvect = np.linspace(0, 1, NZ)
rgrid, phigrid = np.meshgrid(rvect, phivect)

def test_3D_conversion_cst_scale_height():
    #radial powerlaw
    input_data = rgrid**-0.5
    
    output_data = twoD2threeD(arr2d=input_data, scale_height=0.05, zvect=zvect)
    ref = pickle.load(open(test_dir/"ref/2D-3D.p", mode="rb"))
    assert np.all(output_data == ref)
    ##regold with
    #with open(test_dir/"ref/2D-3D.p", mode="wb") as file:
    #    pickle.dump(output_data, file)


def test_3D_conversion_var_scale_height():
    #radial powerlaw
    input_data = rgrid**-0.5
    h = 0.05 * rgrid
    output_data = twoD2threeD(arr2d=input_data, scale_height=h, zvect=zvect)
    ref = pickle.load(open(test_dir/"ref/2D-3D_h.p", mode="rb"))
    assert np.all(output_data == ref)
    #regold with
    #with open(test_dir/"ref/2D-3D_h.p", mode="wb") as file:
    #    pickle.dump(output_data, file)
