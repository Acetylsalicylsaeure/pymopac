import pytest
import sys
from pymopac import MopacInput


def test_without_dependencies(monkeypatch):
    monkeypatch.setitem(sys.modules, 'numpy', None)
    monkeypatch.setitem(sys.modules, 'rdkit', None)
    monkeypatch.setitem(sys.modules, 'matplotlib', None)

    xyz = """C          0.94359       -0.07160       -0.03505
C          2.45567       -0.07160       -0.03505
H          0.55932        0.87684       -0.42205
H          0.55932       -0.21067        0.97983
H          0.55932       -0.88098       -0.66292
H          2.83994        0.73778        0.59283
H          2.83994        0.06746       -1.04993
H          2.83994       -1.02005        0.35196"""

    inp = MopacInput(xyz)
    out = inp.run()
    assert 57.44*0.9 < float(out.COSMO_VOLUME) < 57.44*1.1


def test_correct_failure(monkeypatch):
    monkeypatch.setitem(sys.modules, 'numpy', None)
    monkeypatch.setitem(sys.modules, 'rdkit', None)
    monkeypatch.setitem(sys.modules, 'matplotlib', None)

    with pytest.raises(ImportError) as exc_info:
        inp = MopacInput("CC")
        out = inp.run()
    assert "None in sys" in str(exc_info.value)
