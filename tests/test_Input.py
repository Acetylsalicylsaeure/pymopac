from pymopac import MopacInput
from rdkit import Chem


def test_AddHs():
    inp = MopacInput("CCC", AddHs=True, preopt=True)
    assert len([x for x in inp.mol.GetAtoms()]) == 11
