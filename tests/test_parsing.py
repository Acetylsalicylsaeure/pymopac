import pytest
from pymopac import MopacInput
from rdkit import Chem


def test_modelheader():
    comment = "# test"
    inp = MopacInput("CCC", comment=comment, addHs=True,
                     preopt=True, verbose=True, aux=False)
    print(inp)
    out = inp.run()
    assert out.header == "PM7"
    assert out.comment == comment


def test_modelheader_aux():
    comment = "# test"
    inp = MopacInput("CCC", comment=comment, addHs=True,
                     preopt=True, verbose=True)
    print(inp)
    out = inp.run()
    assert out.header == "PM7 AUX"
    assert out.comment == comment


def test_HF():
    outfile = MopacInput("[H][F]").run()
    assert float(outfile["MOLECULAR_WEIGHT"]) == 20.0063
    assert float(outfile.__dict__["NO._OF_FILLED_LEVELS"]) == 4
    iop = float(outfile["IONIZATION_POTENTIAL"])
    iop_target = 15.795581
    assert iop_target*0.95 < iop < iop_target*1.05
