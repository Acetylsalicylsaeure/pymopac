import rdkit
import pytest
import subprocess
from pymopac import MOPAC_PATH


def test_rdkit():
    v_rdkit = rdkit.__version__
    assert isinstance(v_rdkit, str)


def test_mopac_import():
    assert MOPAC_PATH
