from pymopac.classes import GeometryToMol
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit


def test_smi():
    mol = GeometryToMol("CC")
    assert isinstance(mol, Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2


def test_mol():
    rdmol = Chem.MolFromSmiles("CC")
    mol = GeometryToMol(rdmol)
    assert isinstance(mol, Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2


def test_preoptimized_mol():
    rdmol = Chem.MolFromSmiles("CC")
    AllChem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol)
    AllChem.MMFFOptimizeMolecule(rdmol)
    mol = GeometryToMol(rdmol)
    assert isinstance(mol, rdkit.Chem.rdchem.Mol)
    assert len(mol.GetAtoms()) == 2
