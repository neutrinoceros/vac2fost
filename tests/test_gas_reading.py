import os
import pathlib
import shutil
import subprocess
import pytest

import f90nml
from astropy.io import fits
from vac2fost import main as app
from vac2fost.vac2fost import DETECTED_MCFOST_VERSION

test_dir = pathlib.Path(__file__).absolute().parent
outdir = test_dir / "output/test_read_gas_density"
if outdir.is_dir():
    shutil.rmtree(outdir)
conf_file = test_dir/"sample/vac2fost_conf_quick.nml"
itf = app(conf_file, output_dir=outdir,
          dust_bin_mode="dust-only", read_gas_density=True, mcfost_verbose=True)


@pytest.mark.incremental #each test is run only if the previous one passed
class TestGasReading:
    def test_read_gas_density_header(self):
        header = fits.open(itf.io.OUT.filepath)[0].header
        cards_as_dicts = {c.keyword: c.value for c in header.cards}
        assert ("read_gas_density", 1) in cards_as_dicts.items()

    def test_read_gas_density_shape(self):
        gas_density = fits.open(itf.io.OUT.filepath)[-1].data
        conf = f90nml.read(conf_file)["mcfost_output"]
        target_shape = tuple([conf[k] for k in ("nr","nz", "nphi")])
        assert gas_density.shape == target_shape

    @pytest.mark.skipif(DETECTED_MCFOST_VERSION[2] < "35", reason="latest mcfost versions only")
    def test_read_gas_density_is_valid(self):
        """check that mcfost doesn't crash when passed gas density."""
        os.chdir(itf.io.OUT.directory)
        subprocess.check_call(
            f"mcfost mcfost_conf.para -density_file {itf.io.OUT.filepath.name} -3D",
            shell=True,
        )

