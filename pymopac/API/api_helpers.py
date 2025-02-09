import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from ctypes import c_int, c_double
import numpy as np
from .api_types import MopacSystem
from rdkit.Chem import rdDetermineBonds


def lib_finder():
    conda_prefix = os.getenv("CONDA_PREFIX")
    if conda_prefix:
        suspect = f"{conda_prefix}/lib/libmopac.so"
        if os.path.isfile(suspect):
            return suspect

    try:
        output = subprocess.check_output(["ldconfig", "-p"], text=True)
        for line in output.splitlines():
            if any(lib_name in line for lib_name in ["libmopac.so"]):
                return line.split("=>")[-1].strip()
    except (subprocess.SubprocessError, FileNotFoundError):
        pass

    return None


def mol_to_system(mol, charge=0, spin=0, model=0,
                  preopt=True, add_Hs=True, optimize_geometry=True):
    """Convert an RDKit molecule to a MOPAC system

    Args:
        mol: RDKit molecule object
        preopt: if the mol object should by pre-optimized via MMFF
        assHs: if hydrogens should be substituted at empty spots
        charge: Net molecular charge (default: 0)
        spin: Number of unpaired electrons (default: 0)
        model: Semiempirical model (default: 0/PM7)
        optimize_geometry: Whether atoms can move in calculation (default: True)

    Returns:
        MopacSystem object ready for calculation
    """
    # Ensure molecule has explicit hydrogens and 3D coordinates
    if add_Hs:
        mol = Chem.AddHs(mol)
    if not mol.GetNumConformers():
        AllChem.EmbedMolecule(mol, randomSeed=42)
    if preopt:
        AllChem.MMFFOptimizeMolecule(mol)

    try:
        from pymopac import __mopac_version__
        _version = int(__mopac_version__.replace(".", ""))
        if _version <2303:
            _natom_fix = True
        else:
            _natom_fix=False
    except:
        _natom_fix=True
    # Get basic molecular information
    num_atoms = mol.GetNumAtoms()
    # special condition to add a dummy in case of 4 atoms
    # otherwise, mopac test condition force exits with code 1
    if _natom_fix and num_atoms == 4:
        rwmol = Chem.RWMol(mol)
        rwmol.AddAtom(Chem.Atom(99))
        mol = rwmol.GetMol()
        num_atoms += 1
    conf = mol.GetConformer()

    # Create MOPAC system
    system = MopacSystem()  # Changed from API.MopacSystem()

    # Set basic parameters
    system.natom = num_atoms
    system.natom_move = num_atoms if optimize_geometry else 0
    system.charge = charge
    system.spin = spin
    system.model = model
    system.epsilon = 1.0
    system.nlattice = 0
    system.nlattice_move = 0
    system.pressure = 0.0
    system.tolerance = 1.0
    system.max_time = 3600

    # Create atomic number array
    atomic_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    system.atom = (c_int * num_atoms)(*atomic_nums)

    # Create coordinate array
    coords = []
    for atom_idx in range(num_atoms):
        pos = conf.GetAtomPosition(atom_idx)
        coords.extend([pos.x, pos.y, pos.z])
    system.coord = (c_double * (3*num_atoms))(*coords)

    return system


def system_to_mol(system, sanitize=True, charge=0, cov_factor=1.3):
    """Convert a MOPAC system to an RDKit molecule

    Args:
        system: MopacSystem object
        sanitize: Whether to sanitize the molecule (default: True)
        charge: Molecular charge (default: 0)
        cov_factor: Factor to multiply covalent radii (default: 1.3)

    Returns:
        RDKit Mol object
    """
    # Create empty editable mol object
    mol = Chem.RWMol()

    # Create conformer for 3D coordinates
    conf = Chem.Conformer(system.natom)

    # Get coordinates as numpy array
    coords = np.ctypeslib.as_array(system.coord, shape=(system.natom * 3,))

    # Add atoms and their 3D coordinates
    for i in range(system.natom):
        atom = Chem.Atom(system.atom[i])
        atom_idx = mol.AddAtom(atom)
        conf.SetAtomPosition(
            atom_idx, (coords[i*3], coords[i*3 + 1], coords[i*3 + 2]))

    # Add conformer to molecule
    mol.AddConformer(conf)

    # Let RDKit determine connectivity using van der Waals method
    rdDetermineBonds.DetermineConnectivity(mol,
                                           useHueckel=False,
                                           charge=charge,
                                           covFactor=cov_factor,
                                           useVdw=True)

    if sanitize:
        try:
            mol = mol.GetMol()
            Chem.SanitizeMol(mol)
        except Exception as e:
            print(f"Warning: Sanitization failed: {e}")

    return mol
