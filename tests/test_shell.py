'''tests for installed version, called through shell'''

from pathlib import Path
import subprocess
import pytest

from astropy.io import fits as pyfits

from vac2fost import main as app
from vac2fost import __file__ as v2cfile
root = Path(v2cfile).parent


@pytest.mark.incremental #each test is run only if the previous one passed
class TestShellCalling():
    '''require that vac2fost.py be accessible through your $PATH env variable'''

    def test_command_line_call(self):
        print(f'ROOT : {root}')
        exitcode = subprocess.call(
            str(root /'vac2fost.py ') + str(root.parent / 'tests/sample/vac2fost_conf.nml'),
            shell=True
        )
        assert exitcode == 0

    def test_command_line_call_w_offset(self):
        exitcode = subprocess.call(
            root / 'vac2fost.py ' + str(root.parent / 'test/sample/vac2fost_conf.nml -o 2'),
            shell=True
        )
        assert exitcode == 0


    # def test_install(self):
    #     exitcode = subprocess.call('which vac2fost.py', shell=True)
    #     assert exitcode == 0
