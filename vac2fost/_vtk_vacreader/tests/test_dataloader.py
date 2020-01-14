from pathlib import Path

import numpy as np
import pytest
from vtk_vacreader import VacDataLoader, VacDataSorter

basename = Path(__file__).absolute().parent/ 'data/graindrift'
myshape = (512, 128)


def test_load_unexisting_data():
    vdl = VacDataLoader(
        basename=basename,
        shape=myshape
    )
    with pytest.raises(FileNotFoundError):
        vdl[100]

def test_load_existing_data():
    vdl = VacDataLoader(
        basename=basename,
        shape=myshape
    )

    ntest = 0
    vds = VacDataSorter(
        file_name=f"{basename}{str(ntest).zfill(4)}.vtu",
        shape=myshape
    )
    for key in vds.fields:
        assert np.all(vdl[ntest][key] == vds[key])
