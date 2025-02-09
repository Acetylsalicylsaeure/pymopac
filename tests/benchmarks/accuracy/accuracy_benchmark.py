import pymopac
from pymopac import API
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm


def prepare_mol(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return mol
    except Exception as e:
        print(f"Error preparing molecule {smiles}: {e}")
        return None


def predict_MopacInput(mol):
    if mol is None:
        return None
    try:
        inp = pymopac.MopacInput(mol, preopt=False, addHs=False)
        out = inp.run()
        return out.FINAL_HEAT_OF_FORMATION.number / 4.184
    except Exception as e:
        print(e)
        return None


def predict_API(mol):
    if mol is None:
        return None
    try:
        system = API.mol_to_system(mol, preopt=False, add_Hs=False)
        props, state = API.optimize_geometry(system)
        return props.heat
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    df = pd.read_csv("./pddg_smiles.csv")
    api_results = []
    mopac_results = []

    for smiles in tqdm(df['SMILES']):
        mol = prepare_mol(smiles)
        api_results.append(predict_API(mol))
        mopac_results.append(predict_MopacInput(mol))

    df["API"] = api_results
    df["MopacInput"] = mopac_results
    df.to_csv("result.csv", index=False)
