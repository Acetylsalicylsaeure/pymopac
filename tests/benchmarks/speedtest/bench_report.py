from api_bench import benchmark
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors

# List of common molecules ordered roughly by size/complexity
MOLECULES = [
    "F",
    "C",           # Methane
    "CC(=O)O",     # Acetic acid
    "c1ccccc1",
    "CCc1ccccc1CC",
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
]


def run_benchmarks(molecules, repetitions=10):
    """Run benchmarks for each molecule and combine results"""
    results = []
    for smiles in molecules:
        print(f"Benchmarking {smiles}")
        df = benchmark(smiles, repetitions)
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        df['Molecule'] = smiles
        df['Size'] = Descriptors.ExactMolWt(mol)
        df['Atoms'] = mol.GetNumAtoms()
        results.append(df)

    return pd.concat(results, ignore_index=True)


def plot_results(df):
    """Create two performance plots: absolute times and relative performance"""
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 12), height_ratios=[1, 1])

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # , '#9467bd']

    # First subplot: Absolute performance
    for func, color in zip(sorted(df['Function'].unique()), colors):
        func_data = df[df['Function'] == func]

        # Calculate means and CI for each atom count
        stats = func_data.groupby('Atoms')['Time (ms)'].agg(
            ['mean', 'std', 'count']).reset_index()
        stats['ci'] = 1.96 * stats['std'] / np.sqrt(stats['count'])

        # Add jitter to x-values for scatter points
        x_jittered = func_data['Atoms'].values + \
            np.random.uniform(-0.2, 0.2, size=len(func_data))

        # Plot individual points with jitter
        ax1.scatter(x_jittered, func_data['Time (ms)'].values,
                    color=color, alpha=0.3, s=50)

        # Plot mean line with error bars
        ax1.errorbar(stats['Atoms'], stats['mean'],
                     yerr=stats['ci'],
                     color=color,
                     label=func,
                     marker='o',
                     markersize=8,
                     linewidth=2,
                     linestyle='--',
                     capsize=5,
                     capthick=2,
                     zorder=10)

    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_yscale("log")
    ax1.set_xlabel('Number of Atoms')
    ax1.set_ylabel('Time per Measurement (ms)')
    ax1.set_title('Absolute Performance')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Second subplot: Relative performance
    # Calculate mean times for Classical method for each atom count
    classical_means = df[df['Function'] == 'Classical'].groupby('Atoms')[
        'Time (ms)'].mean()

    # Merge relative times into the DataFrame
    df = df.copy()
    df['Relative Time'] = df.apply(
        lambda row: row['Time (ms)'] / classical_means.loc[row['Atoms']]
        if row['Atoms'] in classical_means else np.nan, axis=1
    )

    for func, color in zip(sorted(df['Function'].unique()), colors):
        func_data = df[df['Function'] == func]

        # Calculate means and CI for relative times
        rel_stats = func_data.groupby('Atoms')['Relative Time'].agg(
            ['mean', 'std', 'count']).reset_index()
        rel_stats['ci'] = 1.96 * rel_stats['std'] / np.sqrt(rel_stats['count'])

        # Add jitter to x-values for scatter points
        x_jittered = func_data['Atoms'].values + \
            np.random.uniform(-0.2, 0.2, size=len(func_data))

        # Plot individual relative measurements
        ax2.scatter(x_jittered, func_data['Relative Time'].values,
                    color=color, alpha=0.3, s=50)

        # Plot mean line with error bars
        ax2.errorbar(rel_stats['Atoms'], rel_stats['mean'],
                     yerr=rel_stats['ci'],
                     color=color,
                     label=func,
                     marker='o',
                     markersize=8,
                     linewidth=2,
                     linestyle='--',
                     capsize=5,
                     capthick=2,
                     zorder=10)

    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_xlabel('Number of Atoms')
    ax2.set_ylabel('Relative Time (compared to Classical)')
    ax2.set_title('Relative Performance (normalized to Classical method)')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Add a bit more space between subplots
    plt.tight_layout()
    plt.savefig('benchmark_plot.svg', bbox_inches='tight', dpi=300)
    plt.close()


def write_report(df):
    """Generate markdown report with results"""
    with open('benchmark.md', 'w') as f:
        f.write('# MOPAC Method Performance Benchmark\n\n')
        f.write('![Benchmark Results](benchmark_plot.svg)\n\n')

        # Overall statistics
        f.write('## Method Performance Summary\n\n')
        summary = df.groupby('Function').agg({
            'Time (ms)': ['count', 'mean', 'std', 'min', 'max']
        }).round(3)
        summary.columns = [
            'Count', 'Mean Time (ms)', 'Std Dev (ms)', 'Min Time (ms)', 'Max Time (ms)']
        f.write(summary.to_markdown() + '\n\n')

        # Results by molecule
        f.write('## Statistics by Molecule\n\n')
        for mol in sorted(df['Molecule'].unique(),
                          key=lambda x: Chem.AddHs(Chem.MolFromSmiles(x)).GetNumAtoms()):
            mol_data = df[df['Molecule'] == mol]
            f.write(f'### {mol} ({mol_data["Atoms"].iloc[0]} atoms)\n\n')

            # Calculate statistics for each method
            stats = mol_data.groupby('Function')['Time (ms)'].agg([
                'count', 'mean', 'std', 'min', 'max'
            ]).round(4)
            stats.columns = ['Measurements', 'Mean', 'Std Dev', 'Min', 'Max']
            f.write(stats.to_markdown() + '\n\n')


if __name__ == "__main__":
    # Run benchmarks and generate report
    results_df = run_benchmarks(MOLECULES)
    plot_results(results_df)
    write_report(results_df)
