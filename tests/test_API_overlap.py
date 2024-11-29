from subprocess import check_call
import pytest
import pymopac
from pymopac import API
from rdkit import Chem
from rdkit.Chem import AllChem


def check_overlap(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    infile = pymopac.MopacInput(mol, preopt=False, addHs=False)
    outfile = infile.run()
    base_mol = outfile.toMol()

    system = API.mol_to_system(mol)
    props, state = API.optimize_geometry(system)
    api_mol = API.system_to_mol(system)

    print("Base mol:")
    print(f"Atoms: {[atom.GetSymbol() for atom in base_mol.GetAtoms()]}")
    print(f"Aromatic: {[a.GetIsAromatic() for a in base_mol.GetAtoms()]}")
    print(f"Bonds: {[(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
          b.GetBondType()) for b in base_mol.GetBonds()]}")

    print("\nAPI mol:")
    print(f"Atoms: {[atom.GetSymbol() for atom in api_mol.GetAtoms()]}")
    print(f"Aromatic: {[a.GetIsAromatic() for a in api_mol.GetAtoms()]}")
    print(f"Bonds: {[(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
          b.GetBondType()) for b in api_mol.GetBonds()]}")

    rmsd = Chem.rdMolAlign.GetBestRMS(base_mol, api_mol)
    print(rmsd)
    assert rmsd < 0.1


def test_benzene_overlap():
    check_overlap("c1ccccc1")


def test_butane_overlap():
    check_overlap("CCCC")


def test_essig_overlap():
    check_overlap("CC(=O)O")
