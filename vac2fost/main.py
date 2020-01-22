"""This main routine creates and uses an interface instance."""
from pathlib import Path

from .info import __version__
from .interfaces import AbstractInterface, VtuFileInterface, DatFileInterface
from .logger import v2flogger as log


def main(
    conf_file: Path,  # python 3.8: positional only
    override: dict = None,
    output_dir: Path = None,
    loglevel: int = 30,
    force_preroll=False,
    hydro_file_type="vtu",
) -> AbstractInterface:
    """Transform a .dat/.vtu datfile into a .fits

    conf_file and overrides are passed down to Interface.__init__()
    loglevel values: 10 debug
                     20 info
                     30 warning
                     40 error
                     50 critical
    """
    log.setLevel(loglevel)
    log.debug(f"start vac2fost {__version__} main loop")
    if hydro_file_type.lower() in ["vtu", ".vtu"]:
        Interface = VtuFileInterface
    else:
        Interface = DatFileInterface
    itf = Interface(conf_file, override=override, output_dir=output_dir)
    try:
        while True:
            log.info(f"current input number: {itf.current_num}\t({itf.iter_frac})")
            try:
                itf.load_input_data()
                itf.preroll_mcfost(force=force_preroll)
                itf.write_output()
            except FileNotFoundError as missing_file:
                log.warning(f"missing file: {missing_file}, resuming iteration.")
            itf.advance_iteration()  # set itf.current_num to next value
    except StopIteration:
        log.debug("end vac2fost")
    return itf  # return the interface object for inspection (tests)
