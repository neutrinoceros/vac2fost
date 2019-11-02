"""This main routine creates and uses an Interface instance."""
from pathlib import Path

from vac2fost.info import __version__
from vac2fost.utils import CYAN, BOLD, get_prompt_size, decorated_centered_message
from vac2fost.interfaces import Interface, VerbatimInterface


def main(config_file: Path, verbose=False, **itf_kwargs):
    """Transform a .vtu datfile into a .fits

    config_file and itf_kwargs are passed down to Interface.__init__()
    """

    print(decorated_centered_message(f"start vac2fost {__version__}"))

    InterfaceType = {True: VerbatimInterface, False: Interface}[verbose]
    itf = InterfaceType(config_file, **itf_kwargs)

    while 1:
        mess1 = f"current input number: {itf.current_num}"
        mess2 = f"({itf.iter_count}/{itf.iter_max})"
        print((BOLD+" "*(get_prompt_size()-len(mess1)-len(mess2))).join([mess1, mess2]))
        print("-"*get_prompt_size())
        try:
            itf.load_input_data()
        except FileNotFoundError as err:
            filepath = Path(str(err)).relative_to(Path.cwd())
            itf.warnings.append(f"file not found: {filepath}")
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
    if itf.warnings:
        print()
        itf.display_warnings()

    print(decorated_centered_message("end vac2fost"))

    # return the Interface object for inspection (tests)
    return itf
