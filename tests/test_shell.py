"""test the shell interface"""

import os
import shutil
from pathlib import Path
from subprocess import check_call, run, CalledProcessError
import pytest

import f90nml

#from vac2fost import __file__ as v2cfile
#root = Path(v2cfile).parent

test_dir = Path(__file__).absolute().parent
OUT = test_dir/"output"
root = test_dir.parent / "vac2fost"

if not OUT.is_dir():
    os.mkdir(OUT)


def test_genconf():
    """check that --genconf outputs a NameList compliant string"""
    output_dir = OUT / "test_genconf"
    if not output_dir.exists(): os.mkdir(output_dir)
    outfile = output_dir / "template_vac2fost.nml"
    with open(outfile, mode="wt") as file:
        run(["python", f"{root}/vac2fost.py", "--genconf"], stdout=file, shell=True, check=True)
    with open(outfile, mode="rt") as file:
        assert f90nml.read(file)

@pytest.mark.incremental #each test is run only if the previous one passed
class TestShellCalling():

    def test_command_line_call(self):
        output_dir = OUT / "test_command_line_call"
        if output_dir.is_dir():
            shutil.rmtree(output_dir)
        check_call([
            "python", f"{root}/vac2fost.py",
            f"{test_dir}/sample/vac2fost_conf.nml", f"--output {output_dir}"
        ])

    def test_command_line_call_w_number(self):
        output_dir = OUT / "test_command_line_call_w_number_1"
        if output_dir.is_dir():
            shutil.rmtree(output_dir)
        run([
            "python",
            f"{root}/vac2fost.py",
            f"{test_dir}/sample/vac2fost_conf_quick.nml",
            f"--output {output_dir}",
            "-n 2"
        ], shell=True, check=True)

    def test_command_line_call_w_number_argonly(self):
        output_dir = OUT / "test_command_line_call_w_number_2"
        if output_dir.is_dir():
            shutil.rmtree(output_dir)
        run([
            "python", f"{root}/vac2fost.py",
            f"{test_dir}/sample/vac2fost_conf_quick_no_number.nml",
            f"--output {output_dir}",
            "-n 2"
        ], shell=True, check=True)

    def test_command_line_call_w_number_argonly_zero(self):
        output_dir = OUT / "test_command_line_call_w_number_3"
        if output_dir.is_dir():
            shutil.rmtree(output_dir)
        run([
            "python", f"{root}/vac2fost.py",
            f"{test_dir}/sample/vac2fost_conf_quick_no_number.nml",
            f"-output {output_dir}",
            "-n 0"
        ], shell=True, check=True)

    def test_command_line_call_wo_number_at_all(self):
        output_dir = OUT / "test_command_line_call_wo_number"
        if output_dir.is_dir():
            shutil.rmtree(output_dir)
        with pytest.raises(CalledProcessError):
            check_call([
                "python", f"{root}/vac2fost.py",
                f"{test_dir}/sample/vac2fost_conf_quick_no_number.nml",
                f"--output {output_dir}",
            ], shell=True)

    def test_command_line_call_w_multiple_numbers(self):
        output_dir = OUT / "test_command_line_w_multiple_numbers"
        if output_dir.is_dir():
            shutil.rmtree(output_dir)
        run([
            "python", f"{root}/vac2fost.py",
            f"{test_dir}/sample/vac2fost_conf_quick_no_number.nml",
            f"--output {output_dir}",
            "-n 0 1 2"
        ], shell=True, check=True)
