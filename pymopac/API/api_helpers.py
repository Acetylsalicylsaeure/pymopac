import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from ctypes import c_int, c_double
import numpy as np
from .api_types import MopacSystem
from rdkit.Chem import rdDetermineBonds
import ctypes.util
import glob
import platform
from typing import Optional

def lib_finder(verbose: bool = False) -> Optional[str]:
    """Find the path to libmopac.so with improved search logic."""
    lib_name = "libmopac.so"
    system = platform.system()

    # 1. Check Conda environments first
    conda_paths = []
    if os.getenv("CONDA_PREFIX"):
        conda_prefix = os.getenv("CONDA_PREFIX")
        conda_paths.extend([
            os.path.join(conda_prefix, "lib", lib_name),
            os.path.join(conda_prefix, "lib64", lib_name),
            os.path.join(conda_prefix, "Library", "bin", lib_name),  # Windows
        ])

    # Also check active Python's site-packages
    try:
        import site
        conda_paths.extend(
            os.path.join(p, lib_name) 
            for p in site.getsitepackages()
        )
    except Exception:
        pass

    for path in conda_paths:
        if os.path.isfile(path):
            if verbose:
                print(f"Found in Conda environment: {path}")
            return path

    # 2. Use ctypes.util.find_library
    ctypes_path = ctypes.util.find_library("mopac")
    if ctypes_path:
        if verbose:
            print(f"Found via ctypes: {ctypes_path}")
        return ctypes_path

    # 3. System-specific searches
    search_paths = []
    if system == "Linux":
        # Check LD_LIBRARY_PATH directories
        ld_paths = os.getenv("LD_LIBRARY_PATH", "").split(":")
        search_paths.extend(ld_paths)
        
        # Common Linux library paths
        search_paths.extend([
            "/usr/local/lib",
            "/usr/lib",
            "/usr/lib/x86_64-linux-gnu",
            "/lib",
            "/lib64"
        ])
        
        # Try ldconfig cache
        try:
            ldconfig = subprocess.run(
                ["ldconfig", "-p"],
                capture_output=True,
                text=True,
                check=True
            )
            for line in ldconfig.stdout.splitlines():
                if lib_name in line:
                    parts = line.strip().split()
                    if len(parts) >= 4 and parts[-2] == "=>":
                        path = parts[-1]
                        if os.path.isfile(path):
                            if verbose:
                                print(f"Found via ldconfig: {path}")
                            return path
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

    elif system == "Darwin":  # macOS
        search_paths.extend([
            "/usr/local/lib",
            "/opt/homebrew/lib",
            "/usr/lib",
            os.path.expanduser("~/lib")
        ])
        
        # Check DYLD_LIBRARY_PATH
        dyld_paths = os.getenv("DYLD_LIBRARY_PATH", "").split(":")
        search_paths.extend(dyld_paths)

    elif system == "Windows":
        search_paths.extend([
            os.path.join(os.getenv("SystemRoot", ""), "System32"),
            os.path.join(os.getenv("PROGRAMFILES", ""), "MOPAC"),
            os.path.join(os.getenv("LOCALAPPDATA", ""), "Programs", "MOPAC")
        ])
        lib_name = "mopac.dll"  # Windows uses different naming

    # 4. Search through all potential paths with version globbing
    for path in search_paths:
        if not os.path.isdir(path):
            continue
            
        # Look for versioned files (e.g., libmopac.so.1.2.3)
        pattern = os.path.join(path, f"{lib_name}*")
        for lib_path in glob.glob(pattern):
            if os.path.isfile(lib_path):
                real_path = os.path.realpath(lib_path)
                if verbose:
                    print(f"Found in system path: {real_path}")
                return real_path

    # 5. Final fallback: try default names with version suffix
    fallback_names = [
        lib_name,
        f"{lib_name}.7",  # Common version suffix pattern
        f"{lib_name}.6",
        f"{lib_name}.1"
    ]
    
    for name in fallback_names:
        path = ctypes.util.find_library(name)
        if path and os.path.isfile(path):
            if verbose:
                print(f"Found via fallback search: {path}")
            return path

    if verbose:
        print("Library not found in any standard locations")
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
