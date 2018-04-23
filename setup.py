from setuptools import setup, find_packages
from pathlib import Path
import os

here = Path(__file__).parent.absolute()
project_name = 'vac2fost'

setup(
    name = project_name,
    version        = __import__(project_name).__version__,
    author         = __import__(project_name).__author__,
    author_email   = __import__(project_name).__contact__,
    url = 'https://gitlab.oca.eu/crobert/vac2fost-project',
    description = 'Interface from AMRVAC data format to MCFOST .fits format.',
    license = 'GNU',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.6'
    ],
    keywords = 'interface data-analysis',
    install_requires = [
        'numpy>=1.13.3',
        'astropy>=2.0.2'
        'f90nml>=1.0.2',
        'amrvac-pywrap>=0.0.2', #same author
    ],
    python_requires = '>=3.6',
    packages = find_packages(),
    package_data = ...
)
