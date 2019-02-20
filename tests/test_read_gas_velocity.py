from vac2fost import main as app
from vac2fost.vac2fost import Interface
import pytest

testconf = "tests/sample/vac2fost_conf.nml"
def test2():
    Interface(testconf, read_gas_velocity=True)

def test1():
    with pytest.raises(NotImplementedError):
        app(testconf, read_gas_velocity=True)

def test3():
    from subprocess import check_call, CalledProcessError
    with pytest.raises(CalledProcessError):
        check_call(["vac2fost.py", f"{testconf}", "--read_gas_velocity"])
