"""Test on grain-size differentiate settling mode"""
from shutil import rmtree

from vac2fost.logger import v2flogger as log
from vac2fost import main as app

from conftest import TEST_DATA_DIR, TEST_ARTIFACTS_DIR

outdir = TEST_ARTIFACTS_DIR / "test_settling"
if outdir.is_dir():
    rmtree(outdir)
conf_file = TEST_DATA_DIR / "vac2fost_conf_quick.nml"


def test_default_setting():
    log.setLevel(10)
    itf = app(conf_file, output_dir=outdir)
    assert not itf._use_settling


def test_run_with_settling():
    log.setLevel(10)
    itf = app(
        conf_file,
        output_dir=outdir,
        override={"flags": dict(settling=True)},
    )
    assert itf._use_settling
    itf.preroll_mcfost()
    itf.get_output_ndarray()
