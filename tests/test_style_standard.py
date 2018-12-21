import subprocess
import pathlib

testloc = pathlib.Path(__file__).parent.resolve()
output_file = testloc/"output/pylint.log" 
mainfile = (testloc/'../vac2fost/vac2fost.py').resolve()

subprocess.call(f"pylint {mainfile} > {output_file}", shell=True)


MINIMAL_MARK = 9.0
def test_style_standard():
    with open(output_file, mode='r') as log:
        lines = log.readlines()
    sumup_line = lines[-2]
    mark = float(sumup_line.split()[6].split('/')[0])
    assert mark >= MINIMAL_MARK

