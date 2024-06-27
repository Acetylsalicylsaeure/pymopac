import pytest
from pymopac import MopacInput
from rdkit import Chem


def test_modelheader():
    comment = "# test"
    inp = MopacInput("CCC", comment=comment, AddHs=True,
                     preopt=True, verbose=True)
    print(inp)
    out = inp.run()
    assert out[list(out.keys())[0]] == "PM7"
    assert out[list(out.keys())[1]] == comment


def test_HF():
    outfile = MopacInput("[H][F]").run()
    assert outfile["MOLECULAR WEIGHT"] == 20.0063
    assert outfile.dic["NO. OF FILLED LEVELS"] == 4
    iop = outfile["IONIZATION POTENTIAL"][0]
    iop_target = 15.795581
    assert iop_target*0.95 < iop < iop_target*1.05
    assert isinstance(outfile.mol, Chem.rdchem.Mol)
