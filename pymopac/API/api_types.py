from ctypes import *
from typing import List

class MopacSystem(Structure):
    """
    Defines the atomistic system and MOPAC job options.

    Attributes:
        natom (int): Number of atoms in the system.
        natom_move (int): Number of atoms allowed to move.
        charge (int): Net charge of the system.
        spin (int): Number of spin excitations.
        model (int): Semiempirical model (PM7=0, PM6-D3H4=1, etc.).
        epsilon (float): Dielectric constant for COSMO solvent.
        atom (POINTER(c_int)): Pointer to atomic numbers array (size `natom`).
        coord (POINTER(c_double)): Pointer to atomic coordinates array (size `3*natom`).
        nlattice (int): Number of lattice vectors.
        nlattice_move (int): Number of moveable lattice vectors.
        pressure (float): External hydrostatic pressure (GPa).
        lattice (POINTER(c_double)): Pointer to lattice vectors array (size `3*nlattice`).
        tolerance (float): Relative numerical tolerance.
        max_time (int): Time limit in seconds.
    """
    _fields_ = [
        ("natom", c_int),
        ("natom_move", c_int),
        ("charge", c_int),
        ("spin", c_int),
        ("model", c_int),
        ("epsilon", c_double),
        ("atom", POINTER(c_int)),
        ("coord", POINTER(c_double)),
        ("nlattice", c_int),
        ("nlattice_move", c_int),
        ("pressure", c_double),
        ("lattice", POINTER(c_double)),
        ("tolerance", c_double),
        ("max_time", c_int)
    ]

class MopacProperties(Structure):
    """
    Stores the calculated ground-state properties and job information.
    
    Attributes:
        heat (float): Heat of formation (kcal/mol).
        dipole (float[3]): Dipole moment vector (Debye).
        charge (POINTER(c_double)): Pointer to atomic partial charges array (size `natom`).
        coord_update (POINTER(c_double)): Pointer to updated coordinates array (size `3*natom_move`).
        coord_deriv (POINTER(c_double)): Pointer to heat gradients array (size `3*natom_move`).
        freq (POINTER(c_double)): Pointer to vibrational frequencies array (size `3*natom_move`).
        disp (POINTER(c_double)): Pointer to displacement vectors array (size `3*natom_move, 3*natom_move`).
        bond_index (POINTER(c_int)): Pointer to CSC format bond matrix indices array (size `natom+1`).
        bond_atom (POINTER(c_int)): Pointer to bonded atoms list (size `bond_index[natom]`).
        bond_order (POINTER(c_double)): Pointer to bond orders array (size `bond_index[natom]`).
        lattice_update (POINTER(c_double)): Pointer to updated lattice vectors array (size `3*nlattice_move`).
        lattice_deriv (POINTER(c_double)): Pointer to lattice gradients array (size `3*nlattice_move`).
        stress (float[6]): Stress tensor in Voigt form.
        nerror (int): Number of error messages.
        error_msg (POINTER(POINTER(c_char))): Pointer to error messages array (size `nerror`).
    """
    _fields_ = [
        ("heat", c_double),
        ("dipole", c_double * 3),
        ("charge", POINTER(c_double)),
        ("coord_update", POINTER(c_double)),
        ("coord_deriv", POINTER(c_double)),
        ("freq", POINTER(c_double)),
        ("disp", POINTER(c_double)),
        ("bond_index", POINTER(c_int)),
        ("bond_atom", POINTER(c_int)),
        ("bond_order", POINTER(c_double)),
        ("lattice_update", POINTER(c_double)),
        ("lattice_deriv", POINTER(c_double)),
        ("stress", c_double * 6),
        ("nerror", c_int),
        ("error_msg", POINTER(POINTER(c_char)))
    ]

class MopacState(Structure):
    """
    Represents an electronic state using standard molecular orbitals.

    Attributes:
        mpack (int): Number of packed matrix elements.
        uhf (int): Unrestricted HF flag (0 = restricted, 1 = unrestricted).
        pa (POINTER(c_double)): Pointer to alpha density matrix array (size `mpack`).
        pb (POINTER(c_double)): Pointer to beta density matrix array (size `mpack`, NULL if `uhf=0`).
    """
    _fields_ = [
        ("mpack", c_int),
        ("uhf", c_int),
        ("pa", POINTER(c_double)),
        ("pb", POINTER(c_double))
    ]

class MozymeState(Structure):
    """
    Represents an electronic state using localized molecular orbitals.

    Attributes:
        numat (int): Number of real atoms.
        nbonds (POINTER(c_int)): Pointer to Lewis bonds per atom array (size `numat`).
        ibonds (POINTER(c_int)): Pointer to Lewis-bonded atoms array (size `9*numat`).
        iorbs (POINTER(c_int)): Pointer to orbitals per atom array (size `numat`).
        noccupied (int): Number of occupied molecular orbitals.
        ncf (POINTER(c_int)): Pointer to atoms in occupied localized molecular orbitals (size `noccupied`).
        nvirtual (int): Number of virtual molecular orbitals.
        nce (POINTER(c_int)): Pointer to atoms in virtual localized molecular orbitals (size `nvirtual`).
        icocc_dim (int): Size of `icocc` array.
        icocc (POINTER(c_int)): Pointer to atom indices in occupied LMOs array (size `icocc_dim`).
        icvir_dim (int): Size of `icvir` array.
        icvir (POINTER(c_int)): Pointer to atom indices in virtual LMOs array (size `icvir_dim`).
        cocc_dim (int): Size of `cocc` array.
        cocc (POINTER(c_double)): Pointer to occupied LMO coefficients array (size `cocc_dim`).
        cvir_dim (int): Size of `cvir` array.
        cvir (POINTER(c_double)): Pointer to virtual LMO coefficients array (size `cvir_dim`).
    """
    _fields_ = [
        ("numat", c_int),
        ("nbonds", POINTER(c_int)),
        ("ibonds", POINTER(c_int)),
        ("iorbs", POINTER(c_int)),
        ("noccupied", c_int),
        ("ncf", POINTER(c_int)),
        ("nvirtual", c_int),
        ("nce", POINTER(c_int)),
        ("icocc_dim", c_int),
        ("icocc", POINTER(c_int)),
        ("icvir_dim", c_int),
        ("icvir", POINTER(c_int)),
        ("cocc_dim", c_int),
        ("cocc", POINTER(c_double)),
        ("cvir_dim", c_int),
        ("cvir", POINTER(c_double))
    ]

class ModelNumber:
    """Enumeration of model numbers for semiempirical methods."""
    PM7 = 0
    PM6_D3H4 = 1
    PM6_ORG = 2
    PM6 = 3
    AM1 = 4
    RM1 = 5

