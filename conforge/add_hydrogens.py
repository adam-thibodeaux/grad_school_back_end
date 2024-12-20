import biotite.structure.io.pdb as biotite_pdb
import biotite.structure.io.mol as biotite_mol
import hydride
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

input_struct = "/Users/adam/Downloads/inputs_for_molec_replac/PAR_TFD_TRIAL_1/paritaprevir_tfd_ektdg_1870.pdb"
output_file_name = "/Users/adam/Downloads/inputs_for_molec_replac/PAR_TFD_TRIAL_1/paritaprevir_tfd_ektdg_1870_with_H.pdb"
reference_bond_struct = "/Users/adam/Downloads/paritaprevir_correct_bonds_no_H.mol"
mol_template = Chem.MolFromMolFile(reference_bond_struct, removeHs=True)
mol = Chem.MolFromPDBFile(input_struct, removeHs=True, sanitize=False, proximityBonding=False)
mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
Chem.MolToMolFile(mol, "./temp.mol")
# print(input_pdb)
        

mol = biotite_mol.MOLFile.read("./temp.mol").get_structure()
mol, _ = hydride.add_hydrogen(mol)
new_mol_file = biotite_mol.MOLFile()
new_mol_file.set_structure(mol)
with open("./temp.mol", "w") as file:
    new_mol_file.write(file)
    file.close()
mol = Chem.MolFromMolFile("./temp.mol", removeHs=False)
Chem.MolToPDBFile(mol, output_file_name)

# 
        
            
#         # res = rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
#         # rdForceFieldHelpers.UFFOptimizeMolecule(mol)
#         # print(res)
#


        # mol = biotite_mol.MOLFile.read("./temp.mol").get_structure()
        
        # mol,_ = hydride.add_hydrogen(mol)
        # temp_pdb = biotite_pdb.PDBFile()
        # temp_pdb.set_structure(mol)
        # with open("./temp.pdb", "w") as file: 
        #     temp_pdb.write(file)
        #     file.close()