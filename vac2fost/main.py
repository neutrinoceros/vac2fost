"""This main routine creates and uses an Interface instance."""
from pathlib import Path

from vac2fost.info import __version__
from vac2fost.interfaces import DatFileInterface, VtuFileInterface
from .logger import v2flogger as log


def main(config_file: Path, loglevel: int = 30, input_data_format="vtu", **itf_kwargs):
    """Transform a .vtu datfile into a .fits

    config_file and itf_kwargs are passed down to Interface.__init__()
    loglevel values: 10 debug
                     20 info
                     30 warning
                     40 error
                     50 critical
    """
    log.setLevel(loglevel)
    log.info(f"start vac2fost {__version__} main loop")
    Interface = {"dat": DatFileInterface, "vtu": VtuFileInterface}[input_data_format]
    itf = Interface(config_file, **itf_kwargs)
    while 1:
        log.info(f"current input number: {itf.current_num}\t({itf.iter_count}/{itf.iter_max})")
        try:
            itf.load_input_data()
        except FileNotFoundError as err:
            filepath = Path(str(err)).relative_to(Path.cwd())
            log.warning(f"missing file: {filepath}, attempting to pursue iteration")
            if itf.iter_count == itf.iter_max:
                break
            continue
        itf.write_mcfost_conf_file()

        if itf_kwargs.get("axisymmetry", False):
            itf.gen_rz_slice()
        else:
            itf.gen_2D_arrays() # todo: rename this method
            itf.gen_3D_arrays() # todo: rename this method
        itf.write_output()

        try:
            filepath = itf.io.OUT.filepath.relative_to(Path.cwd())
        except ValueError:
            filepath = itf.io.OUT.filepath
        try:
            itf.advance_iteration() # set itf.current_num to next value
        except StopIteration:
            break

    log.info("end vac2fost")
    # return the Interface object for inspection (tests)
    return itf
