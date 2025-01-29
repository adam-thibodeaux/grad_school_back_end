import rdkit
from rdkit import Chem

def main(input_pdb_path):

# Read the MOL file
    mol = Chem.MolFromPDBFile(input_pdb_path, removeHs=False)

# Perceive chirality
    Chem.rdmolops.AssignStereochemistryFrom3D(mol)

# Iterate over atoms and check chirality
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            print(f"Atom {atom.GetIdx()}, {atom.GetMonomerInfo().GetName().strip()} is chiral with tag {atom.GetChiralTag()}")

if __name__ == "__main__":
    main("/Users/adam/Downloads/inputs_for_molec_replac/PAR_CONFORGE_TRIAL_1/paritaprevir_alpha_conforge_78.pdb")