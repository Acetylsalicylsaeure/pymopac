from ctypes import *
from .api_types import MopacSystem, MopacState, MopacProperties, MozymeState


def setup_lib(lib):
    """Setup function signatures for the MOPAC library"""

    # MOPAC electronic ground state calculation
    lib.mopac_scf.argtypes = [POINTER(MopacSystem), POINTER(
        MopacState), POINTER(MopacProperties)]
    lib.mopac_scf.restype = None

    # MOPAC geometry relaxation
    lib.mopac_relax.argtypes = [POINTER(MopacSystem), POINTER(
        MopacState), POINTER(MopacProperties)]
    lib.mopac_relax.restype = None

    # MOPAC vibrational calculation
    lib.mopac_vibe.argtypes = [POINTER(MopacSystem), POINTER(
        MopacState), POINTER(MopacProperties)]
    lib.mopac_vibe.restype = None

    # MOZYME electronic ground state calculation
    lib.mozyme_scf.argtypes = [POINTER(MopacSystem), POINTER(
        MozymeState), POINTER(MopacProperties)]
    lib.mozyme_scf.restype = None

    # MOZYME geometry relaxation
    lib.mozyme_relax.argtypes = [POINTER(MopacSystem), POINTER(
        MozymeState), POINTER(MopacProperties)]
    lib.mozyme_relax.restype = None

    # MOZYME vibrational calculation
    lib.mozyme_vibe.argtypes = [POINTER(MopacSystem), POINTER(
        MozymeState), POINTER(MopacProperties)]
    lib.mozyme_vibe.restype = None

    # allocate memory for mopac_state
    lib.create_mopac_state.argtypes = [POINTER(MopacState)]
    lib.create_mopac_state.restype = None

    # allocate memory for mozyme_state
    lib.create_mozyme_state.argtypes = [POINTER(MozymeState)]
    lib.create_mozyme_state.restype = None

    # deallocate memory in mopac_properties
    lib.destroy_mopac_properties.argtypes = [POINTER(MopacProperties)]
    lib.destroy_mopac_properties.restype = None

    # deallocate memory in mopac_state
    lib.destroy_mopac_state.argtypes = [POINTER(MopacState)]
    lib.destroy_mopac_state.restype = None

    # deallocate memory in mozyme_state
    lib.destroy_mozyme_state.argtypes = [POINTER(MozymeState)]
    lib.destroy_mozyme_state.restype = None

    # run MOPAC conventionally from an input file
    lib.run_mopac_from_input.argtypes = [c_char_p]
    lib.run_mopac_from_input.restype = c_int


class MopacLib:
    def __init__(self, lib_path="libmopac.so"):
        self.lib = CDLL(lib_path)  # Load the MOPAC shared library
        setup_lib(self.lib)        # Setup function signatures

    # MOPAC electronic ground state calculation
    def mopac_scf(self, system, state, properties):
        self.lib.mopac_scf(byref(system), byref(state), byref(properties))

    def mopac_relax(self, system, state, properties):  # MOPAC geometry relaxation
        self.lib.mopac_relax(byref(system), byref(state), byref(properties))

    def mopac_vibe(self, system, state, properties):  # MOPAC vibrational calculation
        self.lib.mopac_vibe(byref(system), byref(state), byref(properties))

    # MOZYME electronic ground state calculation
    def mozyme_scf(self, system, state, properties):
        self.lib.mozyme_scf(byref(system), byref(state), byref(properties))

    def mozyme_relax(self, system, state, properties):  # MOZYME geometry relaxation
        self.lib.mozyme_relax(byref(system), byref(state), byref(properties))

    def mozyme_vibe(self, system, state, properties):  # MOZYME vibrational calculation
        self.lib.mozyme_vibe(byref(system), byref(state), byref(properties))

    def create_mopac_state(self, state):  # allocate memory for mopac_state
        self.lib.create_mopac_state(byref(state))

    def create_mozyme_state(self, state):  # allocate memory for mozyme_state
        self.lib.create_mozyme_state(byref(state))

    # deallocate memory in mopac_properties
    def destroy_mopac_properties(self, properties):
        self.lib.destroy_mopac_properties(byref(properties))

    def destroy_mopac_state(self, state):  # deallocate memory in mopac_state
        self.lib.destroy_mopac_state(byref(state))

    def destroy_mozyme_state(self, state):  # deallocate memory in mozyme_state
        self.lib.destroy_mozyme_state(byref(state))

    # run MOPAC conventionally from an input file
    def run_mopac_from_input(self, path):
        return self.lib.run_mopac_from_input(path.encode('utf-8'))


__all__ = ['MopacLib']
