'''test the shell interface'''

import os
from pathlib import Path
import subprocess
import pytest

from astropy.io import fits as pyfits

from vac2fost import main as app
from vac2fost import __file__ as v2cfile

root = Path(v2cfile).parent
test_dir = Path(__file__).absolute().parent

@pytest.mark.incremental #each test is run only if the previous one passed
class TestShellCalling():
    '''require that vac2fost.py be accessible through your $PATH env variable'''

    def test_command_line_call(self):
        output_dir = test_dir / 'output/test_command_line_call/'
        comm = ' '.join([
            str(root / 'vac2fost.py'),
            str(test_dir / 'sample/vac2fost_conf.nml'),
            f'-o {output_dir}'
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode == 0

    def test_command_line_call_w_number(self):
        output_dir = test_dir / 'output/test_command_line_call_w_number_1/'
        comm = ' '.join([
            str(root / 'vac2fost.py'),
            str(test_dir / 'sample/vac2fost_conf_quick.nml'),
            f'-o {output_dir}',
            '-n 2'
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode == 0


    def test_command_line_call_w_number_argonly(self):
        output_dir = test_dir / 'output/test_command_line_call_w_number_2/'
        comm = ' '.join([
            str(root / 'vac2fost.py'),
            str(test_dir / 'sample/vac2fost_conf_quick_no_number.nml'),
            f'-o {output_dir}',
            '-n 2'
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode == 0

    def test_command_line_call_w_number_argonly_zero(self):
        output_dir = test_dir / 'output/test_command_line_call_w_number_3/'
        comm = ' '.join([
            str(root / 'vac2fost.py'),
            str(test_dir / 'sample/vac2fost_conf_quick_no_number.nml'),
            f'-o {output_dir}',
            '-n 0'
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode == 0

    def test_command_line_call_wo_number_at_all(self):
        output_dir = test_dir / 'output/test_command_line_call_wo_number/'
        comm = ' '.join([
            str(root / 'vac2fost.py'),
            str(test_dir / 'sample/vac2fost_conf_quick_no_number.nml'),
            f'-o {output_dir}',
        ])
        exitcode = subprocess.call(comm, shell=True)
        assert exitcode != 0 # Excepted to fail


class TestNarrowCases:
    def test_genconf(self):
        output_dir = test_dir / 'output/test_genconf/'
        comm = ' '.join([
            str(root / 'vac2fost.py'),
            '--genconf',
            f'--output {output_dir}'
        ])
        subprocess.call(comm, shell=True)
        outfile = output_dir / 'template_vac2fost.nml'
        success = False
        assert Path(outfile).exists()
        os.remove(outfile)
