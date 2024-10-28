from pymopac import MopacInput
from rdkit import Chem


def test_addHs():
    inp = MopacInput("CCC", addHs=True, preopt=True)
    assert len(inp.xyz.split("\n")) == 11
