from pymopac.classes import GeometryToMol
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit


def test_smi():
    mol = GeometryToMol("CC", False, False)
    assert isinstance(mol, Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2


def test_mol():
    rdmol = Chem.MolFromSmiles("CC")
    mol = GeometryToMol(rdmol, False, False)
    assert isinstance(mol, Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2


def test_preoptimized_mol():
    rdmol = Chem.MolFromSmiles("CC")
    mol = GeometryToMol(rdmol, True, True)
    assert isinstance(mol, rdkit.Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 8
