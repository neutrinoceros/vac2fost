"""Elementary internal functions and small dataclasses"""
from pathlib import Path
from os.path import expandvars
from shutil import get_terminal_size
from dataclasses import dataclass



# Terminal shenanigans  =============================================================
def shell_path(pin: str) -> Path:
    """Transform <pin> to a Path, expanding included env variables."""
    return Path(expandvars(str(pin)))

def get_prompt_size():
    """size of command line interface messages sized to window. Caps at 80."""
    cols, _ = get_terminal_size()
    return min(cols, 80)



# Dataclasses =======================================================================
@dataclass
class GridShape:
    """Describe number of cells in cylindrical coordinates in a grid"""
    nr: int
    nphi: int
    nz: int = 1

@dataclass
class DataInfo:
    """Hold basic info about input or output data location and shape"""
    directory: Path
    filename: str
    gridshape: GridShape

    @property
    def filepath(self) -> Path:
        """full path"""
        return self.directory / self.filename

    @property
    def filestem(self) -> str:
        """filename without an extension (suffix)"""
        return str(Path(self.filename).stem)

@dataclass
class IOinfo:
    """Hold input and output data information"""
    IN: DataInfo
    OUT: DataInfo



# Decorators ========================================================================
def parameterized(dec):
    """meta decorator, allow definition of decorators with parameters
    source: https://stackoverflow.com/questions/5929107/decorators-with-parameters"""
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)
        return repl
    return layer

@parameterized
def wait_for_ok(func, mess, lenght=get_prompt_size()-7):
    """decorator, sandwich the function execution with '<mess>  ...' & 'ok'"""
    def modfunc(*args, **kwargs):
        print(mess.ljust(lenght), end="... ", flush=True)
        result = func(*args, **kwargs)
        print("ok")
        return result
    return modfunc
