"""Elementary internal functions and small dataclasses"""
from pathlib import Path
from os.path import expandvars
from dataclasses import dataclass


def shell_path(pin: str) -> Path:
    """Transform <pin> to a Path, expanding included env variables."""
    return Path(expandvars(str(pin)))


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
