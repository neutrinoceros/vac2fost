import pathlib
import numpy as np
import pytest
from vac2fost.vac2fost import Interface, MINGRAINSIZE_µ

test_dir = pathlib.Path(__file__).absolute().parent
output_dir = test_dir / "test_dust_binning_mode"

def test_unrecognized_dbm():
    with pytest.raises(KeyError):
        Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                  output_dir=output_dir,
                  dust_bin_mode="not-a-real-dbm-option")

def test_gas_only():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="gas-only")
    assert itf.dust_binning_mode == "gas-only"
    assert itf.grain_micron_sizes == [MINGRAINSIZE_µ]

def test_dust_only():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="dust-only")
    assert itf.dust_binning_mode == "dust-only"
    np.testing.assert_array_equal(itf.grain_micron_sizes, [1e4, 1e3])

def test_mixed():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="mixed")
    assert itf.dust_binning_mode == "mixed"
    np.testing.assert_array_equal(itf.grain_micron_sizes, [MINGRAINSIZE_µ, 1e4, 1e3])

def test_auto_into_mixed():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                    output_dir=output_dir,
                    dust_bin_mode="auto")
    assert itf.dust_binning_mode == "mixed"

def test_auto_into_gas_only():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                    output_dir=output_dir,
                    dust_bin_mode="auto")
    assert itf.dust_binning_mode == "gas-only"

# impossible setups should raise KeyError
def test_mixed_into_KeyError():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                    output_dir=output_dir,
                    dust_bin_mode="mixed")
    with pytest.raises(KeyError):
        gms = itf.grain_micron_sizes

def test_dust_only_into_KeyError():
    itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                    output_dir=output_dir,
                    dust_bin_mode="dust-only")
    with pytest.raises(KeyError):
        gms = itf.grain_micron_sizes
