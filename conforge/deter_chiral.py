import rdkit
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers

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
    # main("/Users/adam/Downloads/inputs_for_molec_replac/PAR_CONFORGE_TRIAL_1/paritaprevir_alpha_conforge_78.pdb")
#     mol_supplier = Chem.SDMolSupplier('./test.sdf', removeHs=False) 
    mol = Chem.MolFromPDBFile("/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_13/sugammadex_subunit_conforge_5_symmetrized_5.pdb", proximityBonding=False, removeHs=False)
    res = rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=100000)
    print(res)
    Chem.MolToPDBFile(mol, "/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_13/sugammadex_subunit_conforge_5_symmetrized_5_opt.pdb")




# # Loop through each molecule in the file
# count = 0
# for mol in mol_supplier:

#     # Perform operations on the molecule (e.g., calculate descriptors, check properties)

#     # print(Chem.MolToSmiles(mol)) 
#     # Chem.MolToMolFile(mol, f"/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_2/sugammadex_sub_unit_conforge_{count}.mol")
#     # count += 1