from rdkit import Chem
from rdkit.Chem import AllChem
from time import time_ns
import os
import subprocess
try:
    from pymopac import MOPAC_PATH
except:
    MOPAC_PATH = "mopac"


class MopacInput():
    """
    Class that sets up the input file for a MOPAC calculation.


    geometry:
        + SMILES (str)

        + rdkit Mol

        + TODO more
    AddHs: bool
        calls AddHs on the mol object
    preopt: bool
        uses MMFF to optimize the mol structure
    model: str
        all supported MOPAC keywords, e.g. AM1, PM6, ...
    comment: str
        second line in the input file
    run: bool
        directly runs the MOPAC calculation with init
    path:
        if False, a dir is created under /tmp/, else under the given str

    """

    def __init__(self, geometry, AddHs: bool = False, preopt: bool = False,
                 model: str = "PM7", comment: str = "#", run: bool = False,
                 path=False):
        self.model = model
        self.comment = comment
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
        """
        concatenates all class features into the MOPAC input file
        """
        header = self.model
        xyz = Chem.MolToXYZBlock(self.mol).split("\n")[1:]
        xyz = "\n".join(xyz)

        inp = "\n".join([header, self.comment, xyz])
        return inp

    def run(self):
        """
        runs MOPAC as a subprocess
        returns MopacOutput class
        """
        with open(f"{self.tmp_dir}/pymopac.mop", "w") as file:
            file.write(self.inpfile())
        os.chdir(self.tmp_dir)
        process = subprocess.run(
            [MOPAC_PATH, "pymopac.mop"],
            capture_output=True)
        if process.returncode == 0:
            print("MOPAC ran sucessfully")
        else:
            print(process.stderr)
        return MopacOutput(out_path=f"{self.tmp_dir}/pymopac.out",
                           stderr=process.stderr, stdout=process.stdout)


class MopacOutput():
    """
    reads the MOPAC .out file at the given out_path and parses datapoints

    standalone runs possible, but calling via MopacInput().run() recommended.
    """

    def __init__(self, out_path: str, stdout=None, stderr=None):
        self.stdout = stdout
        self.stderr = stderr
        if not os.path.isfile(out_path):
            raise Exception("output file not found")

        with open(out_path, "r") as file:
            self.outfile = file.read()


if __name__ == "__main__":
    inp = MopacInput("CCC", AddHs=True, preopt=True)
    print(inp)
    out = inp.run()
    print(out.outfile)
