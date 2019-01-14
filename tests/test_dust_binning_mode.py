import pathlib
import numpy as np
import pytest
from vac2fost.vac2fost import Interface, MINGRAINSIZE_µ

test_dir = pathlib.Path(__file__).absolute().parent
output_dir = test_dir / "test_dust_binning_mode"


def test_gas_only():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="gas-only")
    grain_micron_sizes = itf.grain_micron_sizes
    assert grain_micron_sizes == [MINGRAINSIZE_µ]

def test_dust_only():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="dust-only")
    grain_micron_sizes = itf.grain_micron_sizes
    print(grain_micron_sizes)
    np.testing.assert_array_equal(grain_micron_sizes, [1e4, 1e3])

def test_mixed():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="mixed")
    grain_micron_sizes = itf.grain_micron_sizes
    np.testing.assert_array_equal(grain_micron_sizes, [MINGRAINSIZE_µ, 1e4, 1e3])

def test_auto_into_mixed():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="auto")
    itf.grain_micron_sizes
    assert itf.dust_binning_mode == "mixed"

def test_auto_into_gas_only():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                    output_dir=output_dir,
                    dust_bin_mode="auto")
    itf.grain_micron_sizes
    assert itf.dust_binning_mode == "gas-only"

# impossible setups should raise KeyError
def test_mixed_into_KeyError():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                    output_dir=output_dir,
                    dust_bin_mode="mixed")
    with pytest.raises(KeyError):
        itf.grain_micron_sizes

def test_dust_only_into_KeyError():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                    output_dir=output_dir,
                    dust_bin_mode="dust-only")
    with pytest.raises(KeyError):
        itf.grain_micron_sizes

def test_unknown_dbm():
    with pytest.raises(KeyError):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                        output_dir=output_dir,
                        dust_bin_mode="blblblblblblblb")
