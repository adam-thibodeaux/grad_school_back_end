from rdkit import Chem
import sys
def main():
    mol = Chem.MolFromSmiles(sys.argv[-1])
    smiles_set = set()
    for _ in range(100):
        smiles_str = Chem.MolToSmiles(mol, doRandom=True)
        smiles_set.add(smiles_str)
    print(list(smiles_set))
    sys.stdout.flush()
    return

if __name__ == "__main__":
    main()