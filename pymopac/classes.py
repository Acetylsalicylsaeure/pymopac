from rdkit import Chem
from rdkit.Chem import AllChem
from time import time_ns
import os
import subprocess


class MopacInput():
    """
    Class that sets up the input file for a MOPAC calculation.
    geometry: 
        SMILES (str)
        rdkit Mol
        TODO more
    AddHs: bool
        calls AddHs on the mol object
    preopt: bool
        uses MMFF to optimize the mol structure
    model: str
        all supported MOPAC keywords, e.g. AM1, PM6, ...
    run: bool
        directly runs the MOPAC calculation with init
    path:
        if False, a dir is created under /tmp/, else under the given str

    """

    def __init__(self, geometry, AddHs: bool = False, preopt: bool = False, model: str = "PM7", run: bool = False,
                 path=False):
        if isinstance(geometry, str):
            try:
                self.mol = Chem.MolFromSmiles(geometry)
            except Exception as e:
                raise RuntimeError("Failed to parse input as SMILES", e)
        if isinstance(geometry, Chem.rdchem.Mol):
            self.mol = geometry

        if AddHs:
            self.mol = AllChem.AddHs(self.mol)
        if preopt:
            AllChem.EmbedMolecule(self.mol)
            AllChem.MMFFOptimizeMolecule(self.mol)

        if not path:
            ts = time_ns()
            self.tmp_dir = f"/tmp/pymopac_{ts}"
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = path
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir)

        if run:
            self.run()

    def inpfile(self):
        pass

    def run(self):
        pass


if __name__ == "__main__":
    inp = MopacInput("CCC", AddHs=True, preopt=True)
    print(inp)
    print(inp.mol)
    print(Chem.MolToXYZBlock(inp.mol))
    print("object doc:", inp.__doc__)
    print("class doc:", MopacInput.__doc__)
