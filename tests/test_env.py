# devnote : can't pickle an interfarce object because it contains vtk objects...
import pathlib
import vac2fost
from vac2fost import main as app
from vac2fost.vac2fost import mcfost_major_version, mcfost_minor_version, DETECTED_MCFOST_VERSION
import pytest

testdir = pathlib.Path(__file__).absolute().parent
output_dir = testdir/'output/test_env'

@pytest.mark.skipif(not vac2fost.vac2fost.colorama,
                    reason="Can't test color printing withou colorama")
def test_with_colorama():
    itf = app(
        str(testdir/'sample/vac2fost_conf_quick.nml'),
        output_dir=output_dir
    )
    itf.display_warnings()

def test_without_colorama():
    vac2fost.vac2fost.colorama = None
    itf = app(
        str(testdir/'sample/vac2fost_conf_quick.nml'),
        output_dir=output_dir
    )
    itf.display_warnings()
