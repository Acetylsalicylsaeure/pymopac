"""
PyMOPAC: A Python interface for MOPAC calculations.

This package provides tools and utilities for working with MOPAC,
a semi-empirical quantum chemistry program.
"""

from .helpers import get_mopac
from .classes import MopacInput, MopacOutput

__version__ = "Beta"
__author__ = "Acetylsalicylsaeure"


__all__ = ['get_mopac', "MopacInput", "MopacOutput"]
