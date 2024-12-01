import tempfile
import os
import numpy as np
import pymopac
from pymopac import API
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import subprocess
import timeit
import time
import gc
import random


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

    measurements = {
        'Function': [],
        'Time (ms)': [],
        'Run': [],
        'Order': []  # Track the order of execution
    }

    # Create full experiment design
    all_runs = []
    for run in range(repetitions):
        # For each repetition, create a shuffled list of all functions
        func_names = list(funcs.keys())
        random.shuffle(func_names)
        all_runs.extend([(name, run) for name in func_names])

    # Warmup phase (one run of each function in random order)
    warmup_funcs = list(funcs.items())
    random.shuffle(warmup_funcs)
    for _, func in warmup_funcs:
        func()

    # Execute the randomized experiment design
    for order, (name, run) in enumerate(all_runs):
        gc.collect()
        time.sleep(0.1)

        timer = timeit.Timer(funcs[name])
        wall_time = timer.timeit(number=1) * 1000  # ms

        measurements['Function'].append(name)
        measurements['Time (ms)'].append(wall_time)
        measurements['Run'].append(run)
        measurements['Order'].append(order)

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
