import subprocess
import pathlib

test_dir = pathlib.Path(__file__).parent.resolve()
output_file = test_dir/"output/pylint.log"
mainfile = (test_dir.parent/"vac2fost/vac2fost.py").resolve()
rcfile = test_dir.parent/".pylintrc"
MINIMAL_MARK = 10.0 #perfect score or nothing

def test_style_standard():
    subprocess.check_call(f"pylint {mainfile} --rcfile {rcfile} > {output_file}", shell=True)
    with open(output_file, mode='r') as log:
        lines = log.readlines()
    sumup_line = lines[-2]
    mark = float(sumup_line.split()[6].split('/')[0])
    assert mark >= MINIMAL_MARK
