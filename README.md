# pyMOPAC

initial test build with basic functionality

wraps MOPAC to conveniently interact with the program from within python scripts

## Dependencies
+ MOPAC in $PATH (some search is done for slightly different names, but older versions might not parse)
+ rdkit

+ numpy<2 ([rdkit issue](https://github.com/rdkit/rdkit/issues/7477))

+ matplotlib (optional for plotting)

## Installation

+ via pip

`pip install pymopac`

+ directly from github (nightly)
```bash
git clone Acetylsalicylsaeure/pymopac
cd pymopac
pip install .
```

## Usage

At the core, this module implements two classes, MopacInput and MopacOutput. These represent the MOPAC .mop input and .out output files, respectively. MopacInput takes a molecular geometry in various supported formats like SMILES or rdkit Mol. Furthermore, various keywords are accepted that define the MOPAC Header. The input file can be sent to MOPAC via the .run() method, which returns a MopacOutput object.

This Object Class parses the .out file and aims to dynamically extract calculation results.

### Minimal working example:
```python
import pymopac


infile = pymopac.MopacInput("CC")
outfile = infile.run()
print(outfile.outfile)
```

### Getting calculated properties
```python
import pymopac


outfile = pymopac.MopacInput("c1ccccc1").run()
print(outfile.keys())
print(outfile["IONIZATION POTENTIAL"])
```

### working with Mol objects
The module internally represents the molecule via rdkit Mol objects. They can serve both as inputs and can be accessed as outputs after having their geometry optimized with MOPAC
```python
from rdkit import Chem
from rdkit.Chem import AllChem
import pymopac
from rdkit.Chem import rdDetermineBonds


mmff_mol = Chem.MolFromSmiles("c1ccc1")
mmff_mol = AllChem.AddHs(mmff_mol)
AllChem.EmbedMolecule(mmff_mol)
AllChem.MMFFOptimizeMolecule(mmff_mol)

outfile = pymopac.MopacInput(mmff_mol).run()
mopac_mol = outfile.mol
# since we don't parse bonds at this stage of development, it is necessary to infer them
rdDetermineBonds.DetermineBonds(mopac_mol)

print(AllChem.GetBestRMS(mmff_mol, mopac_mol))
```

### Run feedback
3 different keywords are implemented, which offer feedback to a MOPAC run

+ `verbose=True`
    prints internal messages from the python module to stdout

+ `stream=True`
    streams the MOPAC .out file to stdout

+ `plot=True`
    uses matplotlib to plot the progress in gradient and heat of formation

## Tests

done for Ubuntu 24 LTS and Fedora 40,

Python 3.9-12

## Roadmap
+ add multiline parsing
+ add parsing for further keywords
