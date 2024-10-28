from os.path import isfile
from rdkit.Chem.AllChem import warnings
from .helpers import xyz_identifier, get_mopac, BlockToXyz, checkOverlap
from .output import MopacOutput
import os
from time import time_ns, sleep
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

    def run(self):
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

    Methods:
        + run()
        runs MOPAC as a subprocess
        returns MopacOutput class

        + getInpFile()
        returns the MOPAC input as a string
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

        self.xyz = self.geoToXyz()

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

    def geoToXyz(self):
        """
        Helper function that tries to infer a geometry from any input and
        returns a xyz block
        """
        # check for xyz block
        xyz_status, xyz_block = xyz_identifier(self.geometry)
        if xyz_status:
            return BlockToXyz(xyz_block)

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

        if self.stream or self.plot:
            process = self.verboseRun()
        else:
            process = self.silentRun()

        return MopacOutput(outfile=self.getOutResult(),
                           stderr=process.stderr, stdout=process.stdout,
                           aux=self.getAuxResult())

    def stream_run(self):
        """
        returns both the mopac process and the stream yielding lines
        """

        with open(self.outpath, 'a'):
            os.utime(self.outpath, None)
        while not os.path.isfile(self.outpath):
            sleep(0.01)

        # start tail process first and create stream generator
        tail_process = subprocess.Popen(["tail", "-f", self.outpath],
                                        stdout=subprocess.PIPE,
                                        universal_newlines=True)

        def pure_stream(mopac_proc, tail_proc):
            try:
                for line in tail_proc.stdout:
                    yield line.strip()
                    if mopac_proc.poll() is not None:
                        break
                    if "MOPAC DONE" in line:
                        break
            finally:
                tail_proc.terminate()
                tail_proc.wait()
                mopac_proc.wait()

        # Start MOPAC process after tail is ready
        mopac_process = subprocess.Popen(["mopac", self.inpath])

        # Create stream with both processes
        stream = pure_stream(mopac_process, tail_process)
        return mopac_process, stream

    def verboseRun(self):
        """
        Sets up the basics for a stream run and contains switches for being
        verbose and streaming. If plot is True, it creates a live-updating
        matplotlib figure of gradient values.
        """
        process, lines = self.stream_run()
        line_counter = 0

        if self.plot:
            import matplotlib.pyplot as plt
            from matplotlib.animation import FuncAnimation

            self.grads = []
            self.heats = []
            self.fig, (self.ax, self.ax_heat) = plt.subplots(2, 1)
            self.line, = self.ax.plot([], [])
            self.line_heat, = self.ax_heat.plot([], [])

            self.ax.set_xlabel('Cycle')
            self.ax.set_ylabel('Gradient')
            self.ax_heat.set_xlabel('Cycle')
            self.ax_heat.set_ylabel('Heat')

            self.ax.set_title('Gradient vs. Cycle')
            self.ax_heat.set_title('Heat vs. Cycle')

            def update_plot(frame):
                self.line.set_data(range(len(self.grads)), self.grads)
                self.line_heat.set_data(range(len(self.heats)), self.heats)
                self.ax.relim()
                self.ax.autoscale_view()
                self.ax_heat.relim()
                self.ax_heat.autoscale_view()
                return self.line, self.line_heat

            self.ani = FuncAnimation(self.fig, update_plot,
                                     frames=None, interval=100, blit=True,
                                     cache_frame_data=False, save_count=100)
            plt.tight_layout()
            plt.show(block=False)

        # Process all lines until the generator is exhausted
        for line in lines:
            if self.stream:
                print(line)
            if self.plot and "GRAD.:" in line:
                spl = line.split()
                i = spl.index("GRAD.:")
                i_heat = spl.index("HEAT:")
                self.grads.append(float(spl[i+1]))
                self.heats.append(float(spl[i_heat+1]))
                self.fig.canvas.flush_events()
            line_counter += 1

        # Generator has finished, wait for process to complete
        process.wait()
        if line_counter < 2:
            print("No lines captured, calculations presumably done too fast")

        if self.plot:
            import matplotlib.pyplot as plt
            plt.show()  # Keep the plot open after the function finishes

        return process

    def silentRun(self):
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
            length = len(result_splitter)
            result_splitter += " " + \
                "\n ".join(self.getInpFile().split("\n")[:2])
            try:
                i = out.index(result_splitter) + length
            except:
                warnings.warn("no output splitter found, proceed with caution")
                i = 0
            result = out[i:]
            return result

    def getAuxResult(self):
        if hasattr(self, "auxpath"):
            if os.path.isfile(self.auxpath):
                with open(self.auxpath, "r") as f:
                    return f.read()
        else:
            return None
