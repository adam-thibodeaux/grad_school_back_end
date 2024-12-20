from rmsd_calc import calc_rmsd
import os

input_dir = "./para_conformers_large_trial_3"
compare_struct_dir = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb"
rmsd_list = []
for pdb_file in os.listdir(input_dir):
    rmsd_list.append(calc_rmsd(input_struct_dir=f"{input_dir}/{pdb_file}",compare_struct_dir=compare_struct_dir))
print(sorted(rmsd_list))