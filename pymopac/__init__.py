"""
PyMOPAC: A Python interface for MOPAC calculations.

This package provides tools and utilities for working with MOPAC,
a semi-empirical quantum chemistry program.
"""

from .helpers import get_mopac, get_version_string

__version__ = "1.0"
__author__ = "Acetylsalicylsaeure"
__mopac_path__ = get_mopac()

if __mopac_path__:
    __mopac_version__ = get_version_string(__mopac_path__)


from .input import MopacInput
from .output import MopacOutput

__all__ = ["MopacInput", "MopacOutput", "get_mopac"]
