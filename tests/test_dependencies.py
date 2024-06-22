import rdkit
import pytest
import subprocess


def test_rdkit():
    v_rdkit = rdkit.__version__
    assert isinstance(v_rdkit, str)


def test_mopac():
    v_mopac = subprocess.run(["mopac", "-V"], capture_output=True)
    assert v_mopac.returncode == 0
