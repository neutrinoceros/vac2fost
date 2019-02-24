from vac2fost import main as app, __file__ as v2ffile
from vac2fost.vac2fost import Interface
import pytest
from pathlib import Path


testdir = Path(__file__).parent

testconf = testdir/"sample/vac2fost_conf.nml"
def test1():
    Interface(testconf, read_gas_velocity=True)

def test2():
    with pytest.raises(NotImplementedError):
        app(testconf, read_gas_velocity=True)
