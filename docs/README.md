# pyMOPAC

initial test build with basic functionality

wraps MOPAC to conveniently interact with the program from within python scripts

## Dependencies
+ MOPAC in $PATH (some search is done for slightly different names, but older versions might not parse)
+ rdkit

+ numpy<2 ([rdkit issue](https://github.com/rdkit/rdkit/issues/7477))

## Usage

At the core, this module implements two classes, MopacInput and MopacOutput. These represent the MOPAC .mop input and .out output files, respectively. MopacInput takes a molecular geometry in various supported formats like SMILES or rdkit Mol. Furthermore, various keywords are accepted that define the MOPAC Header. The input file can be sent to MOPAC via the .run() method, which returns a MopacOutput object.

This Object Class parses the .out file and aims to dynamically extract calculation results.

## Tests

done for Ubuntu 24 LTS and Fedora 40,

Python 3.9-12
