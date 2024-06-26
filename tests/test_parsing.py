import pytest
from pymopac import MopacInput


def test_modelheader():
    comment = "# test"
    inp = MopacInput("CCC", comment=comment, AddHs=True,
                     preopt=True, verbose=True)
    print(inp)
    out = inp.run()
    assert out[list(out.keys())[0]] == "PM7"
    assert out[list(out.keys())[1]] == comment
