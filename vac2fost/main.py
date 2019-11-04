"""This main routine creates and uses an interface instance."""
from pathlib import Path

from .info import __version__
from .interfaces import AbstractInterface, DatFileInterface, VtuFileInterface
from .logger import v2flogger as log


def main(conf_file: Path, # python 3.8: positional only
         override: dict = None,
         loglevel: int = 30,
         input_data_format: str = "vtu") -> AbstractInterface:
    """Transform a .vtu datfile into a .fits

    conf_file and overrides are passed down to Interface.__init__()
    loglevel values: 10 debug
                     20 info
                     30 warning
                     40 error
                     50 critical
    """
    log.setLevel(loglevel)
    log.info(f"start vac2fost {__version__} main loop")
    if input_data_format != "vtu":
        raise NotImplementedError
    Interface = {"dat": DatFileInterface, "vtu": VtuFileInterface}[input_data_format] # wip
    itf = Interface(conf_file, override)
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
        itf.preroll_mcfost()
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
    return itf  # return the interface object for inspection (tests)
