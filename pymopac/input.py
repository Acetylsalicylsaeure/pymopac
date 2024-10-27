from .helpers import optional_imports, xyz_identifier, get_mopac, BlockToXyz, checkOverlap
from .output import MopacOutput
import os
from time import time_ns
import subprocess


try:
    from pymopac import get_mopac
    MOPAC_PATH = get_mopac()
    if MOPAC_PATH is None:
        MOPAC_PATH = "mopac"
        print("MOPAC not found in path, falling back to trying 'mopac'. expect issues.")
except Exception as e:
    print(e)
    MOPAC_PATH = "mopac"


class BaseInput:
    """
    Base class that sets up the basic mechanisms of launching a job and getting
    outputs.
    """

    def __init__(self, **kwargs):
        pass


class MopacInput(BaseInput):
    """
    Class that sets up the input file for a MOPAC calculation.


    geometry:
        + pure xyz block

        + SMILES (str)

        + rdkit Mol

    model: str
        all supported MOPAC keywords, e.g. AM1, PM6, ...

    custom_header: str
        custom string that is attached to other header keywords

    comment: str
        second line in the input file

    path:
        if False, a dir is created under /tmp/, else under the given str

    AddHs: bool
        calls AddHs on the mol object

    preopt: bool
        uses MMFF to optimize the mol structure

    verbose: bool
        if True, prints additional class internal statements

    stream: bool
        if True, streams the outfile to stdout

    plot: bool
        if True, plots progress via matplotlib

    """

    def __init__(self, geometry,
                 model: str = "PM7",
                 custom_header: str = "",
                 comment: str = "#",
                 path=False,
                 addHs: bool = True,
                 preopt: bool = True,
                 verbose: bool = False,
                 stream: bool = False,
                 plot: bool = False,
                 aux: bool = True):
        self.geometry = geometry
        self.model = model
        self.custom_header = custom_header
        self.comment = comment
        self.path = path
        self.addHs = addHs
        self.preopt = preopt
        self.verbose = verbose
        self.stream = stream
        self.plot = plot
        self.aux = aux

        if self.plot:
            self.stream = True

        self.xyz = self.GeoToXyz()

        if not path:
            ts = time_ns()
            self.path = f"/tmp/pymopac_{ts}"
            if not os.path.isdir(self.path):
                os.mkdir(self.path)
                if self.verbose:
                    print(f"created path at {self.path}")
        else:
            if not os.path.isdir(self.path):
                os.mkdir(str(self.path))
                if self.verbose:
                    print(f"created path at {self.path}")

    def GeoToXyz(self):
        """
        Helper function that tries to infer a geometry from any input and
        returns a xyz block
        """
        # check for xyz block
        xyz_status, xyz_block = xyz_identifier(self.geometry)
        if xyz_status:
            return xyz_block

        from rdkit import Chem
        from rdkit.Chem import AllChem
        if isinstance(self. geometry, str):
            try:
                mol = Chem.MolFromSmiles(self.geometry)
            except Exception as e:
                raise RuntimeError("Failed to parse input as SMILES", e)
        else:
            mol = self.geometry

        if self.addHs:
            mol = AllChem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
        if self.preopt:
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol)
            # check if there are some overlapping atoms, reoptimize as needed
            overlapping_pairs = checkOverlap(mol)
            if overlapping_pairs is not []:
                # print(overlapping_pairs)
                conf = mol.GetConformer(0)
                for i, j in overlapping_pairs:
                    import numpy as np
                    wiggle = np.random.normal(0, 0.3, 3)
                    conf.SetAtomPosition(i, conf.GetAtomPosition(i) + wiggle)
                    conf.SetAtomPosition(j, conf.GetAtomPosition(j) - wiggle)
                # print(check_overlap(mol))
                # print(Chem.MolToXYZBlock(mol))
                AllChem.MMFFOptimizeMolecule(mol)
        return BlockToXyz(Chem.MolToXYZBlock(mol))

    def getInpFile(self):
        header = self.model + self.custom_header
        if self.aux:
            header += " AUX"
        return "\n".join([header, str(self.comment), "", str(self.xyz)])

    def run(self):
        """
        runs MOPAC as a subprocess
        returns MopacOutput class
        """
        self.inpath = f"{self.path}/pymopac.mop"
        self.outpath = f"{self.path}/pymopac.out"
        if self.aux:
            self.auxpath = f"{self.path}/pymopac.aux"
        with open(self.inpath, "w") as file:
            file.write(self.getInpFile())
        if self.verbose:
            print(f"input file written to {self.inpath}")

        if self.stream:
            process = self.stream_run()
        else:
            process = self.silent_run()

        return MopacOutput(outfile=self.getOutResult(),
                           stderr=process.stderr, stdout=process.stdout,
                           aux=self.getAuxResult())

    def stream_run(self):
        pass

    def verbose_run(self):
        pass

    def silent_run(self):
        """
        just runs MOPAC, no feedback or streaming
        """
        process = subprocess.run(
            [MOPAC_PATH, self.inpath],
            capture_output=True)
        if process.returncode == 0:
            pass
        else:
            raise Exception(process.stderr)
        return process

    def getOutResult(self):
        with open(self.outpath, "r") as f:
            out = f.read()
            result_splitter = "-------------------------------------------------------------------------------\n"
            result_splitter += "\n".join(self.getInpFile().split("\n")[:2])
            try:
                i = out.index(result_splitter) + len(result_splitter)
            except:
                i = 0
            result = out[i:]
            return result

    def getAuxResult(self):
        with open(self.auxpath, "r") as f:
            return f.read()
