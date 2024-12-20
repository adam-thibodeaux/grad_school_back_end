from rdkit import Chem

mol = Chem.MolFromPDBFile("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_frags/paritaprevir_alpha_frag_4.pdb", removeHs=True, proximityBonding=False)

Chem.MolToMolFile(mol,"/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_frags/temp.mol")
mol = Chem.MolFromMolFile("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_frags/temp.mol")
Chem.MolToPDBFile(mol,"/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_frags/paritaprevir_alpha_frag_4_ordered.pdb")
