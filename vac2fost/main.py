"""This main routine creates and uses an Interface instance."""
from pathlib import Path
import logging
import os
from vac2fost.info import __version__
from vac2fost.utils import CYAN, BOLD, get_prompt_size, decorated_centered_message
from vac2fost.interfaces import Interface, VerbatimInterface


def main(config_file: Path, logto=None, loglevel=logging.WARNING, **itf_kwargs):
    """Transform a .vtu datfile into a .fits

    config_file and itf_kwargs are passed down to Interface.__init__()
    """
    logger_kwargs = dict(level=loglevel, format="%(levelname)s - %(message)s")
    if isinstance(logto, (str, Path)):
        logger_kwargs.update(dict(filename=logto, filemode="w"))
    logging.basicConfig(**logger_kwargs)
    ## devnote this isn't supposed to be called more than once per thread.
    # I don't know how to correct for this in cases where the function is called multiple times in a script...


    itf = Interface(config_file, **itf_kwargs)
    logging.info(f"start vac2fost {__version__}")
    while 1:
        logging.info(f"current input number: {itf.current_num}\t({itf.iter_count}/{itf.iter_max})")
        try:
            filename = itf.load_input_data()
            logging.info(f"sucessfully loaded {filename}")
        except FileNotFoundError as err:
            filepath = Path(str(err)).relative_to(Path.cwd())
            logging.warning(f"missing file: {filepath}, attempting to pursue iteration")
            if itf.iter_count == itf.iter_max:
                break
            continue
        mcfost_conffile = itf.write_mcfost_conf_file()
        logging.info(f"successfully wrote {mcfost_conffile}")

        if itf_kwargs.get("axisymmetry", False):
            itf.gen_rz_slice()
        else:
            itf.gen_2D_arrays() # todo: rename this method
            itf.gen_3D_arrays() # todo: rename this method
        output_file = itf.write_output()
        logging.info(f"successfully wrote {output_file}")

        try:
            filepath = itf.io.OUT.filepath.relative_to(Path.cwd())
        except ValueError:
            filepath = itf.io.OUT.filepath
        try:
            itf.advance_iteration() # set itf.current_num to next value
        except StopIteration:
            break
    # wip
    if itf.warnings:
        for w in itf.warnings:
            logging.warning(w)
        #itf.display_warnings()

    logging.info("end vac2fost")
    # return the Interface object for inspection (tests)
    return itf
