import pymopac
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("[K][K]")
AllChem.EmbedMolecule(mol)
conf = mol.GetConformer(0)
for i in [0, 1]:
    conf.SetAtomPosition(i, [0, 0, 0])


def test_overlap():
    infile = pymopac.MopacInput(mol, preopt=False, aux=False)
    out = infile.run()
    assert isinstance(out, pymopac.MopacOutput)
    assert len(out.result) > 2
    assert isinstance(out.result, str)
