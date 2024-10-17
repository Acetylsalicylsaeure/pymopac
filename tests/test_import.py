import pytest
from pymopac.helpers import optional_imports


def test_ImportWarning():
    with pytest.warns(ImportWarning):
        optional_imports("obviously wrong", globals())


def test_functiontypes():
    with pytest.raises(AttributeError):
        optional_imports(123, globals())


def test_strImport():
    optional_imports("import numpy as np", globals())
    arr = np.array([1, 2, 3])
    assert isinstance(arr, np.ndarray)


def test_listImport():
    optional_imports(["import rdkit", "from rdkit import Chem"], globals())
    assert isinstance(rdkit.__version__, str)
    mol = Chem.MolFromSmiles("C")
    assert mol.GetNumAtoms() == 1
