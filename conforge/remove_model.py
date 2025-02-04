import os
import re
import subprocess
import pathlib
import time
from rdkit import Chem
from rdkit.Chem import rdDistGeom, AllChem, rdForceFieldHelpers
from rdkit.Geometry.rdGeometry import Point3D
import hydride
import biotite.structure.io.mol as biotite_mol
import biotite.structure.io.pdb as biotite_pdb
import biotite
import ast
import math

# dir_1 = "/Users/adam/Downloads/inputs_for_molec_replac/MCMM/Par_ensembles/Par_RCE_1-5A_pdb/"
# dir_2 = "/Users/adam/Downloads/inputs_for_molec_replac/MCMM/Par_ensembles/Par_RCE_1-25A_pdb/"
# fixed_ens_dir needs trailing backslash
fixed_ens_dir = "/Users/adam/Downloads/inputs_for_molec_replac/GRAZ_CONFORGE_TRIAL_1/"
input_mtz = "/Users/adam/Downloads/inputs_for_molec_replac/grazoprevir.mtz"
output_dir = "/Users/adam/Downloads/outputs_from_molec_replac/GRAZ_CONFORGE_TRIAL_1/"
frag_ensemble_prefix = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_frags/paritaprevir_alpha_frag_"
model_regex = re.compile("MODEL.*\n")
end_regex = re.compile("ENDMDL\nEND")

# for f in os.listdir(dir_1):
#     with open(dir_1 + f) as file:
#         file_str = file.read()
#         with  open(fixed_ens_dir + f, "w") as write_file:
#             write_file.write(file_str)
#             write_file.close()
#         file.close()


# for f in os.listdir(dir_2):
#     with open(dir_2 + f) as file:
#         file_str = file.read()
#         with  open(fixed_ens_dir + f, "w") as write_file:
#             write_file.write(file_str)
#             write_file.close()
#         file.close()


times = []
def remove_ancillary(file_str):
    file_str = re.sub(model_regex, "", file_str)
    file_str = re.sub(end_regex, "END", file_str)
    return file_str

def search_with_file(input_mtz, input_pdb, out_dir, identity=1, is_frag=False, frag_ensemble_line="", frag_search_line="", frag_hetatm_line="", output_dir=output_dir, should_debug=False, high_resolution_limit=0.1):
   
    # need this to ensure that double bonds are maintainied and right geometry
   
    # Chem.MolToPDBFile(mol, "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha_no_H.pdb")
    # return
    # try:
        # mol_template = Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds_no_H.mol", removeHs=True)
        # print(input_pdb)
        
        # mol = Chem.MolFromPDBFile(input_pdb, removeHs=True, sanitize=False, proximityBonding=False)
        # mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
        # mol = Chem.AddHs(mol, addCoords=True)
            
        # res = rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
        # rdForceFieldHelpers.UFFOptimizeMolecule(mol)
        # print(res)
        # Chem.MolToPDBFile(mol, "./temp.pdb")


        # mol = biotite_pdb.PDBFile.read("./temp.pdb").get_structure(include_bonds=True)[0]
        # mol.charge = 0
        # print('doanflaksdfnlaksdfnlaksdjflaksjfd')
        # mol,_ = hydride.add_hydrogen(mol)
        
        # temp_pdb = biotite_pdb.PDBFile()
        # temp_pdb.set_structure(mol)
        # with open("./temp.pdb", "w") as file: 
        #     temp_pdb.write(file)
        #     file.close()
        # mol = Chem.MolFromPDBFile("./temp.pdb", removeHs=False, sanitize=False, proximityBonding=False)
        # mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
        # Chem.MolToPDBFile(mol,"./temp.pdb")


        
    # except Exception as e:
    #     print(e)
    #     print("failed adding hydrogen to model, skipping")
    #     return
    start = time.time()
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    phaser_str = ""
    if not is_frag:
        phaser_str = f"""/Applications/ccp4-9/bin/phaser << eof
    MODE MR_AUTO
    HKLIN {input_mtz}
    COMPOSITION ATOM C NUM 9
    COMPOSITION ATOM O NUM 6
    COMPOSITION ATOM S NUM 1
    COMPOSITION ATOM H NUM 14
    COMPOSITION ATOM Na NUM 1
    ENSEMBLE par PDBFILE {input_pdb} IDENTITY 0.95
    ENSEMBLE par DISABLE CHECK ON
    ENSEMBLE par HETATM ON
    FORMFACTORS ELECTRON
    PACK SELECT ALL
    ELLG TARGET 100
    PURGE ROT ENABLE OFF
    PURGE TRA ENABLE OFF
    PURGE RNP ENABLE OFF
    SGALTERNATIVE SELECT ALL
    XYZOUT ON ENSEMBLE ON
    XYZOUT ON PACKING ON 
    TOPFILES 5
    KEYWORDS ON
    ZSCORE USE OFF
    ROOT {output_dir}/{out_dir}_out
    SEARCH ENSEMBLE par
    SEARCH METHOD FULL
    eof"""
    proc = subprocess.run([phaser_str], shell=True, capture_output=True)
    # print(proc.stdout)
    end = time.time()

    #  debugging
    if should_debug:
        with open(output_dir + "/" + "stdout.txt", "w") as file:
            file.write(f"run took {end-start} seconds\n\n\n" + proc.stdout.decode())
            file.close()
    
    return end - start
def main(input_mtz=input_mtz, fixed_ens_dir=fixed_ens_dir, output_dir=output_dir):
    print('here')
    similarity_percents = {4: 0.44, 1:0.15, 3:0.23, 2:0.17}
    combos = [[1], [2], [3], [4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4], [1,2,3], [1,3,4], [1,2,4], [2,3,4],[1,2,3,4]]
    # this is for file based searching
    for f in os.listdir(fixed_ens_dir):
        match_digit_regex = re.compile(r"_\d?_?\d?_?\d")
        for match in re.findall(match_digit_regex, f):
            identity = 0
            for frag_num in match.split("_")[1:]:
                identity += similarity_percents.get(frag_num,0)
        if f[0] == ".": continue
        out_dir = pathlib.Path(f).stem
        #TODO: add in time tracking
        if os.path.isdir(output_dir + out_dir + "/"): 
            print(f"PHASER already complete for {f}, skipping...")
            continue
        print(f'starting PHASER for {f}')
        timer = search_with_file(input_mtz, fixed_ens_dir + f, pathlib.Path(f).stem, identity=identity, output_dir=output_dir)
        times.append(timer)
        print(f"Average time of PHASER run for this set: {sum(times)/len(times)} seconds")
        #TODO: add in description grabbing
        print(f"finished PHASER for {f} in {timer} seconds.")
    return
    # #this is for fragment based searching
    # for combo in combos:
    #     str_combo = [str(x) for x in combo]
    #     frag_ensemble_line = ""
    #     frag_search_line=""
    #     frag_hetatm_line = ""
    #     out_dir = "paritaprevir_alpha_"
    #     for string in str_combo:
    #         # print(f"{fixed_ens_dir}{string}.pdb")
    #         # print(similarity_percents.get(int(string),0))
    #         frag_ensemble_line += f"    ENSEMBLE par{string} PDBFILE {frag_ensemble_prefix}{string}.pdb IDENTITY {similarity_percents.get(int(string),0)}\n"
    #         frag_search_line += f"    SEARCH ENSEMBLE par{string}\n"
    #         frag_hetatm_line += f"    ENSEMBLE par{string} HETATM ON\n"
    #         out_dir += "_"
    #         out_dir += string
    #     if os.path.isdir(output_dir+out_dir):
    #         print(f"PHASER already complete for frag_ensemble {str_combo}, skipping...")
    #         continue
    #     print(f'starting PHASER for frag_ensemble {str_combo}')
    #     timer = search_with_file(input_mtz=input_mtz, out_dir=out_dir,input_pdb="",is_frag=True, frag_ensemble_line=frag_ensemble_line, frag_search_line=frag_search_line, frag_hetatm_line=frag_hetatm_line)
    #     times.append(timer)
    #     print(f"Average time of PHASER run for this set: {sum(times)/len(times)} seconds")
    #     print(f"finished PHASER for frag_ensemble {str_combo} in {timer} seconds.")

def main_tfd(input_pdb, out_dir, output_dir, input_mtz_tfd=input_mtz, should_debug=False, high_resolution_limit = 0.1):
    # print(f"starting phaser for {input_pdb.split('/')[-1]}")
    return search_with_file(input_mtz=input_mtz_tfd, input_pdb=input_pdb, out_dir=out_dir, output_dir=output_dir, should_debug=should_debug, high_resolution_limit = high_resolution_limit)
  
if __name__ == "__main__":
    output_dir = "/Users/adam/Downloads/outputs_from_molec_replac/SUG_CONFORGE_TRIAL_9/"
    input_mtz = "/Users/adam/Downloads/inputs_for_molec_replac/sugammadex.mtz"
    input_dir = "/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_15"
    sorted_dir_list = []
    for f in os.listdir(input_dir):
        if f[0] == ".":
            continue
        # if int(f.split("_")[-1].split(".pdb")[0]) > 350:
            # continue
        sorted_dir_list.append(f)
    print(sorted_dir_list)
    # sorted_dir_list = sorted(sorted_dir_list, key= lambda x : int(x.split("_")[-1].split(".pdb")[0]) )
    timers = []
    for i in range(len(sorted_dir_list)):
        f = sorted_dir_list[i]
        input_pdb =f"{input_dir}/{f}"
        phaser_name = f.split('.pdb')[0]
        new_output_dir = output_dir + phaser_name
        if os.path.exists(new_output_dir):
            print(f"PHASER already complete for {phaser_name}, skipping...")
            continue
        
        print(f'starting PHASER for {phaser_name}')
        timer = main_tfd(input_pdb=input_pdb, out_dir=phaser_name, output_dir=new_output_dir, input_mtz_tfd=input_mtz, high_resolution_limit=1.0, should_debug=False)
        timers.append(timer)
        print(f'finished PHASER for {phaser_name} in {timer} seconds')
        print(f"average time for PHASER runs so far: {sum(timers)/len(timers)} seconds")
        print(f"time remaining: {(sum(timers)/len(timers) * (len(sorted_dir_list) - i))/60} minutes")
    # phaser_str = """/Applications/ccp4-9/bin/phaser << eof
    # MODE MR_AUTO
    # HKLIN /Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.mtz
    # ENSEMBLE par PDBFILE /Users/adam/Downloads/inputs_for_molec_replac/PAR_FRAG_EKTDG_TRIAL_1/ROUND_1_2/paritaprevir_core_only_ektdg_5131.pdb IDENTITY 0.25
    # ENSEMBLE par DISABLE CHECK ON
    # ENSEMBLE par HETATM ON
    # FORMFACTORS ELECTRON
    # PACK KEEP HIGH TFZ ON
    # PURGE ROT ENABLE OFF
    # PURGE TRA ENABLE OFF
    # PURGE RNP ENABLE OFF
    # PACK SELECT ALL
    # PACK CUTOFF 1
    # ELLG TARGET 100
    # SGALTERNATIVE SELECT ALL
    # XYZOUT ON ENSEMBLE ON
    # XYZOUT ON PACKING ON 
    # TOPFILES 5
    # KEYWORDS ON
    # ZSCORE USE OFF
    # ROOT /Users/adam/Downloads/outputs_from_molec_replac/test/test_out
    # SEARCH PRUNE OFF
    # SEARCH ENSEMBLE par
    # eof"""
    # proc = subprocess.run([phaser_str], shell=True, capture_output=True)
    # with open("/Users/adam/Downloads/outputs_from_molec_replac/test/stdout.txt", "w") as file:
    #     file.write(proc.stdout.decode())
    #     file.close()
    # main()
    # tfd_val_path = "/Users/adam/Downloads/inputs_for_molec_replac/PAR_TFD_TRIAL_4_with_H/tfd_out.txt"
    # with open(tfd_val_path, "r") as file:
    #     tfd_vals = ast.literal_eval(file.read())
    #     file.close()
    # input_prefix = "/Users/adam/Downloads/inputs_for_molec_replac/PAR_TFD_TRIAL_4_with_H/paritaprevir_tfd_ektdg_"
    # output_prefix = "/Users/adam/Downloads/outputs_from_molec_replac/PAR_TFD_TRIAL_4_with_H/"
    # for i in range(0,len(tfd_vals)):
    #     tfd_val = tfd_vals[i]
    #     # print(tfd_val)
    #     input_pdb = input_prefix + str(tfd_val[0])
    #     input_pdb += ".pdb"
    #     output_dir = output_prefix + str(tfd_val[0])
    #     if not os.path.isdir(output_dir):
    #         os.mkdir(output_dir)
    #     else:
    #         print(f"phaser already complete for {input_pdb}, skipping...")
    #         continue
    #     print(f"starting phaser for {input_pdb}")
    #     timer = main_tfd(input_pdb=input_pdb, output_dir=output_dir, out_dir=f"paritaprevir_tfd_ekdtg_{tfd_val[0]}")
    #     print(f"phaser completed for {input_pdb} in {timer} seconds")
