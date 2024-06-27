from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
from time import time_ns
import os
import time
import subprocess

try:
    from pymopac import MOPAC_PATH
except Exception as e:
    print(e)
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
    elif isinstance(geometry, rdkit.Chem.rdchem.Mol):
        mol = geometry
    return mol


def ResultFromOutput(outfile: str) -> str:
    break_str = " -------------------------------------------------------------------------------"
    out_l = outfile.split("\n")
    break_index = out_l.index(break_str)+1
    result_l = out_l[break_index:]
    return result_l


def ParseResult(result: str) -> dict:
    """
    takes the result section of the MOPAC outfile, returns targeted parts
    mostly in the format key: (value, unit)
    """
    dic = dict()
    dic["header"] = result[0][1:]
    dic["comment"] = result[1][1:]

    for line in result:
        parsed = ParseLine(line)
        if parsed is not None:
            key, value, unit, store = parsed
            if unit is not None:
                unit = " ".join(unit)
                try:
                    unit = float(unit)
                except Exception as e:
                    pass
                dic[key] = (value, unit)
            else:
                dic[key] = value

    return dic


def ParseLine(line):
    targets = {"FINAL HEAT OF FORMATION =": (-2, None),
               "COSMO AREA              =": (-3, None),
               "COSMO VOLUME            =": (-3, None),
               "GRADIENT NORM           =": (-3, None),
               "IONIZATION POTENTIAL    =": (-2, None),
               "HOMO LUMO ENERGIES (EV) =": (-2, None),
               "NO. OF FILLED LEVELS    =": -1,
               "MOLECULAR WEIGHT        =": 3
               }
    for key in targets.keys():
        if key in line:
            spli = line.split()
            target_key = targets[key]
            if isinstance(target_key, tuple):
                unit_i = target_key[1]
                return (key.strip("=").strip(), float(spli[target_key[0]]),
                        spli[target_key[0]+1:unit_i], None)
            elif isinstance(target_key, int):
                return (key.strip("=").strip(), float(spli[target_key]), None,
                        None)
    return None


def ExtractMol(result: str):
    start_i = None
    for i, n in enumerate(result):
        if "CARTESIAN COORDINATES" in n:
            start_i = i
    if start_i is None:
        return Exception("cartesian coordinates not found")
    start_i += 2
    end_i = result[start_i:].index("")
    end_i += start_i
    XYZRaw = result[start_i:end_i]
    XYZBlock = str(len(XYZRaw)) + "\n\n"

    for line in XYZRaw:
        XYZBlock += " ".join(line.split()[1:]) + "\n"

    mol = Chem.MolFromXYZBlock(XYZBlock)

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

    custom_header: str
        custom string that is attached to other header keywords

    comment: str
        second line in the input file

    path:
        if False, a dir is created under /tmp/, else under the given str

    verbose: bool
        if True, prints additional class internal statements

    stream: bool
        if True, streams the outfile to stdout

    plot: bool
        if True, plots progress via matplotlib

    """

    def __init__(self, geometry,
                 AddHs: bool = True,
                 preopt: bool = True,
                 model: str = "PM7",
                 custom_header: str = "",
                 comment: str = "#",
                 path=False,
                 verbose: bool = False,
                 stream: bool = False,
                 plot: bool = False):
        self.model = model
        self.custom_header = custom_header
        self.comment = comment
        self.verbose = verbose
        self.stream = stream
        self.plot = plot

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
        header = self.model + self.custom_header
        xyz = Chem.MolToXYZBlock(self.mol).split("\n")[1:]
        xyz = "\n".join(xyz)

        inp = "\n".join([header, self.comment, xyz])
        return inp

    def run(self, verbose: bool = None,
            stream: bool = None,
            plot: bool = None):
        """
        runs MOPAC as a subprocess
        returns MopacOutput class

        verbose: bool
            if True, prints additional class internal statements

        stream: bool
            if True, streams the outfile to stdout

        plot: bool
            if True, plots progress via matplotlib
        """
        if verbose is None:
            verbose = self.verbose
        if stream is None:
            stream = self.stream
        if plot is None:
            plot = self.plot

        infile = f"{self.tmp_dir}/pymopac.mop"
        outfile = f"{self.tmp_dir}/pymopac.out"
        with open(infile, "w") as file:
            file.write(self.inpfile())
        if verbose:
            print(f"input file written to {infile}")

        if not verbose and not stream and not plot:
            process = self.silent_run(infile=infile)
        else:
            process = self.verbose_run(infile=infile, verbose=verbose,
                                       stream=stream, plot=plot)

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
                    stream: bool = False,
                    plot: bool = False):
        """
        Sets up the basics for a stream run and contains switches for being
        verbose and streaming. If plot is True, it creates a live-updating
        matplotlib figure of gradient values.
        """
        process, lines = self.stream_run(infile)
        line_counter = 0

        if plot:
            import matplotlib.pyplot as plt
            from matplotlib.animation import FuncAnimation

            self.grads = []
            self.heats = []
            self.fig, (self.ax, self.ax_heat) = plt.subplots(
                2, 1)
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

        for line in lines:
            if stream:
                print(line)
            if plot:
                if "GRAD.:" in line:
                    spl = line.split()
                    i = spl.index("GRAD.:")
                    i_heat = spl.index("HEAT:")
                    self.grads.append(float(spl[i+1]))
                    self.heats.append(float(spl[i_heat+1]))
                    # self.fig.canvas.draw_idle()
                    self.fig.canvas.flush_events()
            line_counter += 1

        if line_counter < 2:
            print("No lines captured, calculations presumably done too fast")

        if plot:
            plt.show()  # Keep the plot open after the function finishes

        return process

    def verbose_run_dep(self, infile: str, verbose: bool = False,
                        stream: bool = False,
                        plot: bool = False):
        """
        sets up the basics for a stream run and contains switches for being
        verbose and streaming
        """
        process, lines = self.stream_run(infile)
        line_counter = 0
        if plot:
            self.grads = []

        for line in lines:
            if stream:
                print(line)
            if plot:
                if "GRAD.:" in line:
                    spl = line.split()
                    i = spl.index("GRAD.:")
                    self.grads.append(spl[i+1])
                    self.fig
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

    outputs are available raw under self.outfile or parsed as a dictionary
    mostly in the format key: (value, unit) under self.dic

    some dictionary functionality is directly available for the class, i.e.
    keys can directly queried via self[key]

    standalone runs possible, but calling via MopacInput().run() recommended.
    """

    def __init__(self, out_path: str, stdout=None, stderr=None):
        self.stdout = stdout
        self.stderr = stderr
        if not os.path.isfile(out_path):
            raise Exception("output file not found")

        with open(out_path, "r") as file:
            self.outfile = file.read()

        self.result = ResultFromOutput(self.outfile)
        self.dic = ParseResult(self.result)
        self.mol = ExtractMol(self.result)

    def keys(self):
        return self.dic.keys()

    def values(self):
        return self.dic.values()

    def items(self):
        return self.dic.items()

    def __getitem__(self, key):
        return self.dic[key]


if __name__ == "__main__":
    inp = MopacInput("[H][F]", verbose=True)
    print(inp)
    out = inp.run()
    for item in out.items():
        print(item)
