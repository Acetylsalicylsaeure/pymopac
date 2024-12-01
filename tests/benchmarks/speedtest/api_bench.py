import tempfile
import os
import numpy as np
import pymopac
from pymopac import API
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import cProfile
import pstats
import subprocess


def benchmark(smiles, repetitions):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    inpfile = pymopac.MopacInput(mol, preopt=False, addHs=False)
    inpstring = inpfile.getInpFile()
    funcs = {
        'MopacInput': lambda: bench_MopacInput(mol),
        'MopacInput_disk': lambda: bench_MopacInput_disk(mol),
        'API': lambda: bench_API(mol),
        'Classical': lambda: bench_classical(inpstring)
    }

    # Create list to store individual measurements
    measurements = {
        'Function': [],
        'Time (ms)': [],
        'Run': [],
        'Function Calls': []
    }

    for name, func in funcs.items():
        for run in range(repetitions):
            profiler = cProfile.Profile()
            profiler.enable()
            func()
            profiler.disable()
            stats_obj = pstats.Stats(profiler)

            # Store individual measurement
            measurements['Function'].append(name)
            measurements['Time (ms)'].append(stats_obj.total_tt * 1000)
            measurements['Run'].append(run)
            measurements['Function Calls'].append(stats_obj.total_calls)

    # Cleanup
    if os.path.exists('./tmp'):
        for f in os.listdir('./tmp'):
            os.remove(os.path.join('./tmp', f))

    return pd.DataFrame(measurements).round(6)


def create_unique_path():
    import uuid
    return os.path.join('./tmp', f'temp_{uuid.uuid4()}.mop')


"""
def bench_API_fromFile(inpstring):
    path = create_unique_path()
    with open(path, "w") as f:
        f.write(inpstring)
    API.run_from_file(path)
"""


def bench_API_fromFile(inpstring):
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp.write(inpstring.encode())
        temp_path = temp.name
    API.run_from_file(temp_path)
    # os.remove(temp_path)


def bench_classical(inpstring):
    path = create_unique_path()
    with open(path, "w") as f:
        f.write(inpstring)
    subprocess.run(['mopac', path], capture_output=True)


def bench_MopacInput(mol):
    inp = pymopac.MopacInput(mol, preopt=False, addHs=False, aux=False)
    out = inp.run()


def bench_MopacInput_disk(mol):
    inp = pymopac.MopacInput(mol, preopt=False, addHs=False, aux=False,
                             path="./tmp")
    out = inp.run()


def bench_API(mol):
    system = API.mol_to_system(mol)
    props, state = API.optimize_geometry(system)


if __name__ == "__main__":
    result = benchmark("C", 5)
    print("\nSample of individual measurements:")
    print(result.head(10))

    print("\nSummary statistics:")
    summary = result.groupby('Function').agg({
        'Time (ms)': ['count', 'mean', 'std', 'min', 'max']
    }).round(3)
    print(summary)
