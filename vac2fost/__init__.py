"""A conversion facility for MPI-AMRVAC (.vtu) to MCFOST (.fits)"""

from vac2fost.info import __version__, __author__, __contact__
from vac2fost.interfaces import VtuFileInterface, DatFileInterface
from vac2fost.main import main
