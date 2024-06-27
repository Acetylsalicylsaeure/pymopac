import rdkit
import pytest
import subprocess
from pymopac.helpers import get_mopac


def test_rdkit():
    v_rdkit = rdkit.__version__
    assert isinstance(v_rdkit, str)


def test_mopac_import():
    assert get_mopac()
