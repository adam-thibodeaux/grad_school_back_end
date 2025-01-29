# from rmsd_calc import calc_rmsd
# import os

# input_dir = "./para_conformers_large_trial_3"
# compare_struct_dir = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb"
# rmsd_list = []
# for pdb_file in os.listdir(input_dir):
#     rmsd_list.append(calc_rmsd(input_struct_dir=f"{input_dir}/{pdb_file}",compare_struct_dir=compare_struct_dir))
# print(sorted(rmsd_list))
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.rdForceFieldHelpers as force_field
print(f"adding hydrogens to refined structure")

mol_ref = Chem.MolFromPDBFile("/Users/adam/Downloads/outputs_from_molec_replac/PAR_CUSTOM_CONF_TRIAL_2_20A/ROUND_10/paritaprevir_torsion_angle_perturb_28_out.1_refine_001.pdb", removeHs=False)
to_remove_atoms = ["H13","H14"]
double_bond_atoms = ["C9","C10"]
fixed_dihedrals = [["H13","C9","C10","H14"]]
fixed_angles = [["C8","C9","C10"],["C11","C10","C9"]]
double_bond_atoms_idx = []
constrained_part = Chem.MolFromPDBFile("/Users/adam/Downloads/outputs_from_molec_replac/PAR_CUSTOM_CONF_TRIAL_2_20A/ROUND_10/paritaprevir_torsion_angle_perturb_28_out.1_refine_001.pdb", removeHs=False)
constrained_part = Chem.EditableMol(constrained_part)
atom_idx_dict = {}

for atom_to_remove in to_remove_atoms:
    for atom in constrained_part.GetMol().GetAtoms():
        if atom.GetMonomerInfo().GetName().strip() == atom_to_remove:
            constrained_part.RemoveAtom(atom.GetIdx())
            break




# AllChem.AssignBondOrdersFromTemplate(Chem.MolFromMolFile("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_fixed_bonds.mol"),mol_ref)
# constrained_part = Chem.DeleteSubstructs(constrained_part,Chem.MolFromSmiles("C=C"))
Chem.MolToPDBFile(constrained_part.GetMol(),"/Users/adam/Downloads/outputs_from_molec_replac/PAR_CUSTOM_CONF_TRIAL_2_20A/ROUND_10/paritaprevir_torsion_angle_perturb_28_out.1_refine_001_testing.pdb")
coord_dict = {}
conf = mol_ref.GetConformer()
atom_idx_dict = {}
for atom in mol_ref.GetAtoms():
    atom_name = atom.GetMonomerInfo().GetName().strip()
    atom_idx = atom.GetIdx()
    if atom_name in double_bond_atoms:
        double_bond_atoms_idx.append(atom_idx)
    if atom_name not in to_remove_atoms:
        coord_dict[atom.GetIdx()] = conf.GetAtomPosition(atom_idx)
    atom_idx_dict[atom_name] = atom_idx
bond = mol_ref.GetBondBetweenAtoms(double_bond_atoms_idx[0],double_bond_atoms_idx[1])
bond.SetBondType(Chem.BondType.DOUBLE)

mp = force_field.MMFFGetMoleculeProperties(mol_ref)
ff = force_field.MMFFGetMoleculeForceField(mol_ref, mp)
for coord in coord_dict.items():
    ff.MMFFAddPositionConstraint(coord[0], 0, 1.e4)
for fixed_dihedral in fixed_dihedrals:
    ff.MMFFAddTorsionConstraint(atom_idx_dict[fixed_dihedral[0]],atom_idx_dict[fixed_dihedral[1]],atom_idx_dict[fixed_dihedral[2]],atom_idx_dict[fixed_dihedral[3]],False, 0.0,0.0,1.e4)
# for angle in fixed_angles:
    # angle_idx_1 = atom_idx_dict[angle[0]]
    # angle_idx_2 = atom_idx_dict[angle[1]]
    # angle_idx_3 = atom_idx_dict[angle[2]]
    # ff.MMFFAddAngleConstraint(angle_idx_1,angle_idx_2,angle_idx_3, False, 120.0,120.0,1.e4)
ff.Minimize(maxIts=10000000)
# conf_id = rdDistGeom.EmbedMolecule(mol_ref, ignoreSmoothingFailures=True, coordMap=coord_dict, maxAttempts=10,)
# print(conf_id)
Chem.MolToPDBFile(mol_ref,"/Users/adam/Downloads/outputs_from_molec_replac/PAR_CUSTOM_CONF_TRIAL_2_20A/ROUND_10/paritaprevir_torsion_angle_perturb_28_out.1_refine_001_testing_2.pdb")