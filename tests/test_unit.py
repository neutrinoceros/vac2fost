import pathlib
import pickle
import pytest
import numpy as np
from astropy.io import fits as pyfits

from vac2fost.vac2fost import twoD2threeD, Interface, MCFOSTUtils
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
    np.testing.assert_array_equal(output_data, ref)
    ##regold with
    #with open(test_dir/"ref/2D-3D.p", mode="wb") as file:
    #    pickle.dump(output_data, file)

def test_3D_conversion_var_scale_height():
    #radial powerlaw
    input_data = rgrid**-0.5
    h = 0.05 * rgrid
    output_data = twoD2threeD(arr2d=input_data, scale_height=h, zvect=zvect)
    ref = pickle.load(open(test_dir/"ref/2D-3D_h.p", mode="rb"))
    np.testing.assert_array_equal(output_data, ref)
    #regold with
    #with open(test_dir/"ref/2D-3D_h.p", mode="wb") as file:
    #    pickle.dump(output_data, file)

def test_3D_conversion_large_grid():
    NR, NPHI, NZ = 23, 61, 7
    rvect = np.linspace(1, 10, NR)
    phivect = np.linspace(0, 2*np.pi, NPHI)
    zvect = np.linspace(0, 1, NZ)
    rgrid, phigrid = np.meshgrid(rvect, phivect)

    #radial powerlaw
    input_data = rgrid**-0.5
    output_data = twoD2threeD(arr2d=input_data, scale_height=0.01, zvect=zvect)
    ref = pickle.load(open(test_dir/"ref/2D-3D_large.p", mode="rb"))
    np.testing.assert_array_equal(output_data, ref)
    ##regold with
    #with open(test_dir/"ref/2D-3D_large.p", mode="wb") as file:
    #    pickle.dump(output_data, file)

def test_3D_conversion_real_usecase():
    ref = pickle.load(open(test_dir/"ref/main_out.p", mode="rb"))
    h = 0.01 * ref['output_grid']['rg']
    zvect = np.linspace(0, 5, 25)
    for a2D, a3D in zip(ref["new_2D_arrays"], ref["new_3D_arrays"]):
        new_3D_arr = twoD2threeD(arr2d=a2D, scale_height=h, zvect=zvect)
        np.testing.assert_allclose(new_3D_arr, a3D, rtol=1e-15)

def test_path_reading():
   """Check that AMRVAC config file can be correctly assessed with a relative "origin" argument"""
   Interface(config_file=test_dir/'sample/vac2fost_conf_nonaxisym.nml', dbg=True)
