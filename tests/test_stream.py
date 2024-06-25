import pymopac
import pytest


def test_streamend(capsys):
    target = "**          Digital Object Identifier (DOI): 10.5281/zenodo.6511958          **"
    target = "ENDED NORMALLY"
    infile = pymopac.MopacInput(
        "CCCCCCC", AddHs=True, stream=True)
    outfile = infile.run()

    captured = capsys.readouterr()
    print(captured.out)
    assert target in captured.out


def test_short_stream(capsys):
    target = "No lines captured, calculations presumably done too fast"

    infile = pymopac.MopacInput("C", AddHs=True, stream=True)
    infile.run()

    captured = capsys.readouterr()
    assert target in captured.out
