import os
from pathlib import Path

from vtk_vacreader import VacDataSorter

datafile = Path(__file__).absolute().parent / 'data/hd142527_dusty0029.vtu'
datafile = str(datafile)

myshape = (512, 128)

def test_load():
    dh = VacDataSorter(datafile)
    assert 1

def test_dict():
    dh = VacDataSorter(datafile)
    dh['rho']
    dh['m1']

def test_resahping():
    dh = VacDataSorter(datafile, shape=myshape)
    assert dh['rho'].shape == myshape
