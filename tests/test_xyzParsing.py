import pytest
from pymopac import MopacInput


def test_simple_xyzblock():
    string = """8

C          0.99799        0.00723       -0.09956
C          2.51006        0.00723       -0.09956
H          0.61372       -0.85259        0.45726
H          0.61371       -0.04508       -1.12259
H          0.61371        0.91935        0.36665
H          2.89433       -0.90490       -0.56577
H          2.89433        0.05954        0.92347
H          2.89433        0.86704       -0.65638"""
    inp = MopacInput(string)
    runfile = inp.getInpFile()
    target_inp = MopacInput("CC")
    target_runfile = target_inp.getInpFile()
    assert len(runfile.split("\n")) == len(target_runfile.split("\n"))
