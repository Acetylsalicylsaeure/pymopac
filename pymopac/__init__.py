"""
PyMOPAC: A Python interface for MOPAC calculations.

This package provides tools and utilities for working with MOPAC,
a semi-empirical quantum chemistry program.
"""

from .helpers import get_mopac
import subprocess

__version__ = "1.0"
__author__ = "Acetylsalicylsaeure"
__mopac_path__ = get_mopac()

if __mopac_path__:
    version_string = subprocess.run(
        [__mopac_path__, "--version"], capture_output=True, text=True)
    version_string = str(version_string.stdout.strip())
    __mopac_version__ = version_string.split()[2]


from .input import MopacInput
from .output import MopacOutput

__all__ = ["MopacInput", "MopacOutput", "get_mopac"]
