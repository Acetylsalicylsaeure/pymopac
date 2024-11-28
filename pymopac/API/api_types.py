from ctypes import *
from typing import List


class MopacSystem(Structure):
    """Data that defines the atomistic system and MOPAC job options"""
    _fields_ = [
        ("natom", c_int),              # number of atoms
        ("natom_move", c_int),         # number of atoms allowed to move
        ("charge", c_int),             # net charge
        ("spin", c_int),               # number of spin excitations
        ("model", c_int),              # semiempirical model (PM7=0, PM6-D3H4=1, etc)
        ("epsilon", c_double),         # dielectric constant for COSMO solvent
        ("atom", POINTER(c_int)),      # atomic numbers [natom]
        ("coord", POINTER(c_double)),  # coordinates [3*natom]
        ("nlattice", c_int),          # number of lattice vectors
        ("nlattice_move", c_int),     # number of moveable lattice vectors
        ("pressure", c_double),        # external hydrostatic pressure (GPa)
        ("lattice", POINTER(c_double)),  # lattice vectors [3*nlattice]
        ("tolerance", c_double),       # relative numerical tolerances
        ("max_time", c_int)           # time limit in seconds
    ]


class MopacProperties(Structure):
    """Calculated ground-state properties and job info"""
    _fields_ = [
        ("heat", c_double),               # heat of formation (kcal/mol)
        ("dipole", c_double * 3),         # dipole moment vector (Debye)
        ("charge", POINTER(c_double)),    # atomic partial charges [natom]
        # updated coordinates [3*natom_move]
        ("coord_update", POINTER(c_double)),
        ("coord_deriv", POINTER(c_double)),  # heat gradients [3*natom_move]
        # vibrational frequencies [3*natom_move]
        ("freq", POINTER(c_double)),
        # displacement vectors [3*natom_move,3*natom_move]
        ("disp", POINTER(c_double)),
        # CSC format bond matrix indices [natom+1]
        ("bond_index", POINTER(c_int)),
        # bonded atoms list [bond_index[natom]]
        ("bond_atom", POINTER(c_int)),
        ("bond_order", POINTER(c_double)),  # bond orders [bond_index[natom]]
        # updated lattice vectors [3*nlattice_move]
        ("lattice_update", POINTER(c_double)),
        # lattice gradients [3*nlattice_move]
        ("lattice_deriv", POINTER(c_double)),
        ("stress", c_double * 6),         # stress tensor in Voigt form
        ("nerror", c_int),                # number of error messages
        ("error_msg", POINTER(POINTER(c_char)))  # error messages [nerror][*]
    ]


class MopacState(Structure):
    """Electronic state using standard molecular orbitals"""
    _fields_ = [
        ("mpack", c_int),              # number of packed matrix elements
        # unrestricted HF flag (0=restricted, 1=unrestricted)
        ("uhf", c_int),
        ("pa", POINTER(c_double)),     # alpha density matrix [mpack]
        # beta density matrix [mpack], NULL if uhf=0
        ("pb", POINTER(c_double))
    ]


class MozymeState(Structure):
    """Electronic state using localized molecular orbitals"""
    _fields_ = [
        ("numat", c_int),              # number of real atoms
        ("nbonds", POINTER(c_int)),    # Lewis bonds per atom [numat]
        ("ibonds", POINTER(c_int)),    # Lewis-bonded atoms [9,numat]
        ("iorbs", POINTER(c_int)),     # orbitals per atom [numat]
        ("noccupied", c_int),          # number of occupied MOs
        ("ncf", POINTER(c_int)),       # atoms in occupied LMOs [noccupied]
        ("nvirtual", c_int),           # number of virtual MOs
        ("nce", POINTER(c_int)),       # atoms in virtual LMOs [nvirtual]
        ("icocc_dim", c_int),          # size of icocc array
        # atom indices in occupied LMOs [icocc_dim]
        ("icocc", POINTER(c_int)),
        ("icvir_dim", c_int),          # size of icvir array
        # atom indices in virtual LMOs [icvir_dim]
        ("icvir", POINTER(c_int)),
        ("cocc_dim", c_int),           # size of cocc array
        ("cocc", POINTER(c_double)),   # occupied LMO coefficients [cocc_dim]
        ("cvir_dim", c_int),           # size of cvir array
        ("cvir", POINTER(c_double))    # virtual LMO coefficients [cvir_dim]
    ]

# Convenience enums for models


class ModelNumber:
    PM7 = 0
    PM6_D3H4 = 1
    PM6_ORG = 2
    PM6 = 3
    AM1 = 4
    RM1 = 5
