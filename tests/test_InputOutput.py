import pymopac
import pytest


def test_WaterPipeline():
    inp = pymopac.MopacInput(geometry="O", AddHs=True)
    outp = inp.run()
    outfile = outp.outfile
    print(outfile)
    assert isinstance(outfile, str)
    assert "ended normally" in outfile.lower()


def test_WaterPipeline_noHs():
    inp = pymopac.MopacInput(geometry="O", AddHs=False)
    outp = inp.run()
    outfile = outp.outfile
    print(outfile)
    assert isinstance(outfile, str)
    assert "ended normally" in outfile.lower()
