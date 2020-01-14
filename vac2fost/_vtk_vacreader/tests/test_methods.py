from pathlib import Path
from vtk_vacreader import VacDataSorter

def test_iter():
    datafile = Path(__file__).resolve().parent / 'data/graindrift0000.vtu'
    data = VacDataSorter(str(datafile), (512,128))
    for (i,v1), (j,v2) in zip(data.fields.items(), data):
        assert i == j
        assert (v1 == v2).all()
    
