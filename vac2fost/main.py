"""This main routine creates and uses an Interface instance."""
from pathlib import Path

import os
from vac2fost.info import __version__
from vac2fost.utils import CYAN, BOLD, get_prompt_size, decorated_centered_message
from vac2fost.interfaces import Interface, VerbatimInterface
from .logger import log as log


def main(config_file: Path, loglevel: int = logging.WARNING, **itf_kwargs):
    """Transform a .vtu datfile into a .fits

    config_file and itf_kwargs are passed down to Interface.__init__()
    loglevel values: 10 debug
                     20 info
                     30 warning
                     40 error
                     50 critical
    """
    log.setLevel(loglevel)
    itf = Interface(config_file, **itf_kwargs)

    log.info(f"start vac2fost {__version__}")
    while 1:
        log.info(f"current input number: {itf.current_num}\t({itf.iter_count}/{itf.iter_max})")
        try:
            filename = itf.load_input_data()
            log.info(f"successfully loaded {filename}")
        except FileNotFoundError as err:
            filepath = Path(str(err)).relative_to(Path.cwd())
            log.warning(f"missing file: {filepath}, attempting to pursue iteration")
            if itf.iter_count == itf.iter_max:
                break
            continue
        mcfost_conffile = itf.write_mcfost_conf_file()
        log.info(f"successfully wrote {mcfost_conffile}")

        if itf_kwargs.get("axisymmetry", False):
            itf.gen_rz_slice()
        else:
            itf.gen_2D_arrays() # todo: rename this method
            itf.gen_3D_arrays() # todo: rename this method
        output_file = itf.write_output()
        log.info(f"successfully wrote {output_file}")

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
