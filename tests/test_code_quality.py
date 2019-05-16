import subprocess
import pathlib

test_dir = pathlib.Path(__file__).parent.resolve()
output_file = test_dir/"output/pylint.log"
install_dir = test_dir.parent
rcfile = install_dir/".pylintrc"
MINIMAL_MARK = 10.0 #perfect score or nothing

def test_style_standard():
    subprocess.check_call(f"pylint {install_dir}/vac2fost/*py --rcfile {rcfile} > {output_file}", shell=True)
    with open(output_file, mode='r') as log:
        lines = log.readlines()
    sumup_line = lines[-2]
    mark = float(sumup_line.split()[6].split('/')[0])
    assert mark >= MINIMAL_MARK
