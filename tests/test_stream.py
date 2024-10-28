import pymopac
import pytest


def test_streamend(capsys):
    # target = "**          Digital Object Identifier (DOI): 10.5281/zenodo.6511958          **"
    target = "ENDED NORMALLY"
    infile = pymopac.MopacInput(
        "CCCCCCCC", addHs=True, stream=True)
    outfile = infile.run()

    captured = capsys.readouterr()
    print(captured.out)
    assert target in captured.out
