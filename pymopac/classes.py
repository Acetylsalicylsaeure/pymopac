from rdkit import Chem
from rdkit.Chem import AllChem
from time import time_ns
import os
import time
import subprocess
import threading

try:
    from pymopac import MOPAC_PATH
except:
    MOPAC_PATH = "mopac"


def GeometryToMol(geometry):
    """
    Helper function that tries to infer a geometry from any input and returns a
    rdkti Mol
    """
    if isinstance(geometry, str):
        try:
            mol = Chem.MolFromSmiles(geometry)
        except Exception as e:
            raise RuntimeError("Failed to parse input as SMILES", e)
    elif isinstance(geometry, Chem.rdchem.Mol):
        mol = geometry
    return mol


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

    path:
        if False, a dir is created under /tmp/, else under the given str

    verbose: bool
        if True, prints additional class internal statements

    stream: bool
        if True, streams the outfile to stdout

    """

    def __init__(self, geometry, AddHs: bool = False, preopt: bool = False,
                 model: str = "PM7", comment: str = "#", path=False,
                 verbose: bool = False, stream: bool = False):
        self.model = model
        self.comment = comment
        self.verbose = verbose
        self.stream = stream

        self.mol = GeometryToMol(geometry)

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
                if self.verbose:
                    print(f"created path at {self.tmp_dir}")
        else:
            self.tmp_dir = path
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir)

    def inpfile(self):
        """
        concatenates all class features into the MOPAC input file
        """
        header = self.model
        xyz = Chem.MolToXYZBlock(self.mol).split("\n")[1:]
        xyz = "\n".join(xyz)

        inp = "\n".join([header, self.comment, xyz])
        return inp

    def run(self, verbose: bool = None, stream: bool = None):
        """
        runs MOPAC as a subprocess
        returns MopacOutput class

        verbose: bool
            if True, prints additional class internal statements

        stream: bool
            if True, streams the outfile to stdout
        """
        if verbose is None:
            verbose = self.verbose
        if stream is None:
            stream = self.stream

        infile = f"{self.tmp_dir}/pymopac.mop"
        outfile = f"{self.tmp_dir}/pymopac.out"
        with open(infile, "w") as file:
            file.write(self.inpfile())
        if verbose:
            print(f"input file written to {infile}")

        if not verbose and not stream:
            process = self.silent_run(infile=infile)
        else:
            process = self.verbose_run(infile=infile, verbose=verbose,
                                       stream=stream)

        return MopacOutput(out_path=outfile,
                           stderr=process.stderr, stdout=process.stdout)

    def silent_run(self, infile: str):
        """
        just runs MOPAC, no feedback or streaming
        """
        process = subprocess.run(
            [MOPAC_PATH, infile],
            capture_output=True)
        if process.returncode == 0:
            pass
        else:
            raise Exception(process.stderr)
        return process

    def verbose_run(self, infile: str, verbose: bool = False,
                    stream: bool = False):
        """
        sets up the basics for a stream run and contains switches for being
        verbose and streaming
        """
        process, lines = self.stream_run(infile)
        if stream is True:
            line_counter = 0
            for line in lines:
                print(line)
                line_counter += 1
            if line_counter < 2:
                print(
                    "no lines captured, calculations presumably done too fast")
        return process

    def stream_run(self, infile):
        """
        returns both the mopac process and the stream yielding lines
        """
        outfile = os.path.splitext(infile)[0] + ".out"
        # Start the MOPAC process
        mopac_process = subprocess.Popen(["mopac", infile])

        # Wait for the output file to be created
        while not os.path.exists(outfile):
            time.sleep(0.1)

        def pure_stream():
            # Start the tail process to stream the output file
            tail_process = subprocess.Popen(["tail", "-f", outfile],
                                            stdout=subprocess.PIPE,
                                            universal_newlines=True)
            try:
                # Stream the output
                for line in tail_process.stdout:
                    yield line.strip()

                    # Check if MOPAC process has finished
                    if mopac_process.poll() is not None:
                        break
                    if "MOPAC DONE" in line:
                        break
            finally:
                # Ensure we terminate the tail process
                tail_process.terminate()
                tail_process.wait()
                # Wait for MOPAC to finish if it hasn't already
                mopac_process.wait()
        return mopac_process, pure_stream()


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
    inp = MopacInput("CCC", AddHs=True, preopt=True, verbose=True)
    print(inp)
    out = inp.run()
