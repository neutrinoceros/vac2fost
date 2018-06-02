'''tests for installed version, called through shell'''

from pathlib import Path
import subprocess
import pytest

from astropy.io import fits as pyfits

from vac2fost import main as app
from vac2fost import __file__ as v2cfile

root = Path(v2cfile).parent
here = Path(__file__).absolute().parent

@pytest.mark.incremental #each test is run only if the previous one passed
class TestShellCalling():
    '''require that vac2fost.py be accessible through your $PATH env variable'''

    def test_command_line_call(self):
        output_dir = here / 'output/test_command_line_call/'
        comm = ' '.join([
            str(root /'vac2fost.py'),
            str(here / 'sample/vac2fost_conf.nml'),
            f'-d {output_dir}'
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode == 0

    def test_command_line_call_w_offset(self):
        output_dir = here / 'output/test_command_line_call_w_offset/'
        comm = ' '.join([
            str(root /'vac2fost.py'),
            str(here / 'sample/vac2fost_conf.nml'),
            f'-d {output_dir}',
            '-o 2'
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode == 0
