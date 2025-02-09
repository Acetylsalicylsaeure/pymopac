from .api_helpers import lib_finder, mol_to_system, system_to_mol
from .api_types import MopacSystem, MopacProperties, MopacState, MozymeState
from .api_subroutines import MopacLib


__lib_path__ = lib_finder()
_lib = MopacLib(__lib_path__)

# Create clean interface functions


def calculate_scf(system: MopacSystem, state: MopacState = None, properties: MopacProperties = None):
    """Perform MOPAC electronic ground state calculation"""
    if state is None:
        state = MopacState()
        _lib.create_mopac_state(state)
    if properties is None:
        properties = MopacProperties()

    _lib.mopac_scf(system, state, properties)
    return properties, state


def optimize_geometry(system: MopacSystem, state: MopacState = None, properties: MopacProperties = None):
    """Perform MOPAC geometry optimization"""
    if state is None:
        state = MopacState()
        _lib.create_mopac_state(state)
    if properties is None:
        properties = MopacProperties()

    _lib.mopac_relax(system, state, properties)
    return properties, state


def calculate_frequencies(system: MopacSystem, state: MopacState = None, properties: MopacProperties = None):
    """Perform MOPAC vibrational analysis"""
    if state is None:
        state = MopacState()
        _lib.create_mopac_state(state)
    if properties is None:
        properties = MopacProperties()

    _lib.mopac_vibe(system, state, properties)
    return properties, state

# Add MOZYME equivalents


def calculate_mozyme_scf(system: MopacSystem, state: MozymeState = None, properties: MopacProperties = None):
    """Perform MOZYME electronic ground state calculation"""
    if state is None:
        state = MozymeState()
        _lib.create_mozyme_state(state)
    if properties is None:
        properties = MopacProperties()

    _lib.mozyme_scf(system, state, properties)
    return properties, state


def optimize_mozyme_geometry(system: MopacSystem, state: MozymeState = None, properties: MopacProperties = None):
    """Perform MOZYME geometry optimization"""
    if state is None:
        state = MozymeState()
        _lib.create_mozyme_state(state)
    if properties is None:
        properties = MopacProperties()

    _lib.mozyme_relax(system, state, properties)
    return properties, state


def calculate_mozyme_frequencies(system: MopacSystem, state: MozymeState = None, properties: MopacProperties = None):
    """Perform MOZYME vibrational analysis"""
    if state is None:
        state = MozymeState()
        _lib.create_mozyme_state(state)
    if properties is None:
        properties = MopacProperties()

    _lib.mozyme_vibe(system, state, properties)
    return properties, state


def run_from_file(filepath: str) -> int:
    """Run MOPAC calculation from input file"""
    return _lib.run_mopac_from_input(filepath)


# Expose classes and functions
__all__ = [
    'MopacSystem',
    'MopacProperties',
    'MopacState',
    'MozymeState',
    'calculate_scf',
    'optimize_geometry',
    'calculate_frequencies',
    'calculate_mozyme_scf',
    'optimize_mozyme_geometry',
    'calculate_mozyme_frequencies',
    'run_from_file',
    'mol_to_system',
    'system_to_mol'
]

if __name__ == "__main__":
    print(__lib_path__)
