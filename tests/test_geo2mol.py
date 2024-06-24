from pymopac.classes import GeometryToMol
import pytest
from rdkit import Chem


def test_smi():
    mol = GeometryToMol("CC")
    assert isinstance(mol, Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2


def test_mol():
    rdmol = Chem.MolFromSmiles("CC")
    mol = GeometryToMol(rdmol)
    assert isinstance(mol, Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2
