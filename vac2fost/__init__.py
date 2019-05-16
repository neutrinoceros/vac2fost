"""A conversion facility for MPI-AMRVAC (.vtu) to MCFOST (.fits)

Disclaimer:
  This package is using Python3.7 syntax/features and will not be made backward
  compatible with older versions of Python.
"""

from vac2fost.info import __version__, __author__, __contact__
from vac2fost.interfaces import Interface, VerbatimInterface
from vac2fost.main import main
