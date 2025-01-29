from rdkit.Chem import rdForceFieldHelpers, AllChem
import rdkit.Chem as Chem

mol = Chem.MolFromPDBFile("/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CUSTOM_CONF_TRIAL_2/ROUND_1/PHASER/paritaprevir_perturb_custom_dihed_11/paritaprevir_perturb_custom_dihed_11_out.1.pdb", removeHs=False, proximityBonding=False)
mol_template = Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds.mol", removeHs=False)
mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=10000)
Chem.MolToPDBFile(mol,"/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CUSTOM_CONF_TRIAL_2/ROUND_1/PHASER/paritaprevir_perturb_custom_dihed_11/paritaprevir_perturb_custom_dihed_11_out.1_optimized.pdb")


