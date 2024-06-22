"""
PyMOPAC: A Python interface for MOPAC calculations.

This package provides tools and utilities for working with MOPAC,
a semi-empirical quantum chemistry program.
"""

from .helpers import get_mopac

__version__ = "early_alpha"
__author__ = "Acetylsalicylsaeure"


# Call get_mopac to find the MOPAC binary
MOPAC_PATH = get_mopac()

if MOPAC_PATH:
    print(f"MOPAC binary found at: {MOPAC_PATH}")
else:
    print("Warning: MOPAC binary not found in PATH")


__all__ = ['get_mopac', 'MOPAC_PATH']
