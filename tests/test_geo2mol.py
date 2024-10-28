from pymopac import MopacInput
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit


def test_smi():
    inp = MopacInput("CC", addHs=False)
    assert isinstance(inp.xyz, str)
    assert len(inp.xyz.split("\n")) == 2


def test_mol():
    rdmol = Chem.MolFromSmiles("CC")
    inp = MopacInput(rdmol, addHs=False)
    assert isinstance(inp.xyz, str)
    assert len(inp.xyz.split("\n")) == 2


def test_preoptimized_mol():
    rdmol = Chem.MolFromSmiles("CC")
    inp = MopacInput(rdmol, addHs=False)
    assert isinstance(inp.xyz, str)
    assert len(inp.xyz.split("\n")) == 2
