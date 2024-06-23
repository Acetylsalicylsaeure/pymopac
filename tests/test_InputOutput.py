import pymopac
import pytest


def test_WaterPipeline():
    inp = pymopac.MopacInput(geometry="O", AddHs=True, preopt=True)
    outp = inp.run()
    outfile = outp.outfile
    print(outfile)
    assert isinstance(outfile, str)
    assert "ended normally" in outfile.lower()
