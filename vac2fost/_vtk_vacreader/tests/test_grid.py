import os
from pathlib import Path

import numpy as np
from vtk_vacreader import VacDataSorter

datafile = Path(__file__).resolve().parent / 'data/graindrift0000.vtu'

datafile = str(datafile)
myshape = (512, 128)

def test_grid():
    dh = VacDataSorter(datafile, myshape)
    mstar = 2.2
    vphi_d = (dh['m2d1']/dh['rhod1']).mean(axis=1)
    vk = np.sqrt(mstar / dh.get_ticks(axis=0))
    assert abs((vk-vphi_d)/vphi_d).max() < 2e-5
