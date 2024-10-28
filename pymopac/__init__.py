"""
PyMOPAC: A Python interface for MOPAC calculations.

This package provides tools and utilities for working with MOPAC,
a semi-empirical quantum chemistry program.
"""

from .helpers import get_mopac
from .input import MopacInput
from .output import MopacOutput

__version__ = "Stable"
__author__ = "Acetylsalicylsaeure"


__all__ = ["MopacInput", "MopacOutput", "get_mopac"]
