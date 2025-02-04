import os
from io import StringIO
import pandas as pd
import re
import rmsd_calc
import matplotlib.pyplot as plt
input_dir = "/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CUSTOM_CONF_TRIAL_4_20A/"
reference_struct_path = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta_reordered.pdb"
dir_dict = {}
for f in os.listdir(input_dir):
    if "OLD" in f or "ROUND" not in f: continue
    for inner_f in os.listdir(input_dir + f):
        if "PHASER" not in inner_f: continue
        for inner_inner_f in os.listdir(input_dir + f + "/" + inner_f):
            if "summary" not in inner_inner_f: continue
            with open(input_dir + f + "/" + inner_f + "/" + inner_inner_f) as file:
                file_clean = re.sub("P\s+[0-9]*\s+[0-9]*\s+[0-9]*", "", file.read())
            
                file_clean = re.sub("space_group", "", file_clean)
                df = pd.read_csv(StringIO(file_clean), sep="\s+", on_bad_lines="skip")
                # df.sort_values("tfz", inplace=True, ascending=False)
                df.sort_values(["llg","tfz"], inplace=True, ascending=[False,False])
                dir_dict[f] = df
                file.close()


prefix = "ROUND_"
list_dir_1 = []
list_dir_2 = []
list_dir_3 = []
for i in range(len(dir_dict.values())):
    print(f"THIS IS FOR ROUND {i}!!")
    key = prefix + str(i)
    print(dir_dict[key])
    data_1 = dir_dict[key].iloc[0]
    data_2 = dir_dict[key].iloc[1]
    data_3 = dir_dict[key].iloc[2]
    names = [x["name"].split("_out")[0] for x in [data_1]]
    count = 1
    while data_2["name"].split("_out")[0] in names and count < len(dir_dict[key]):
        print(data_2)
        data_2 = dir_dict[key].iloc[count]
        count += 1
    names.append(data_2["name"].split("_out")[0])
    while data_3["name"].split("_out")[0] in names and count < len(dir_dict[key]):
        data_3 = dir_dict[key].iloc[count]
        count += 1
    names.append(data_3["name"].split("_out")[0])
    print(names)
    try:
        rmsd_1 = rmsd_calc.calc_rmsd(input_dir + key + "/" + data_1["name"].split("_out")[0] + ".pdb", reference_struct_path)
        rmsd_2 = rmsd_calc.calc_rmsd(input_dir + key + "/" + data_2["name"].split("_out")[0] + ".pdb", reference_struct_path)
        rmsd_3 = rmsd_calc.calc_rmsd(input_dir + key + "/" + data_3["name"].split("_out")[0] + ".pdb", reference_struct_path)
        top_1 = {"top solution 1": rmsd_1, "top solution 2": rmsd_2, "top solution 3": rmsd_3}
    except Exception as e:
        print(e)
        # top_1 = {"top solution 1": 0, "top solution 2": 0, "top solution 3": 0}
        # try:
        #     rmsd_1 = rmsd_calc.calc_rmsd_pdb(input_dir + key + "/" + data_1["name"].split("_out")[0] + ".pdb", reference_struct_path)
        #     rmsd_2 = rmsd_calc.calc_rmsd_pdb(input_dir + key + "/" + data_2["name"].split("_out")[0] + ".pdb", reference_struct_path)
        #     rmsd_3 = rmsd_calc.calc_rmsd_pdb(input_dir + key + "/" + data_3["name"].split("_out")[0] + ".pdb", reference_struct_path)
        #     top_1 = {"top solution 1": rmsd_1, "top solution 2": rmsd_2, "top solution 3": rmsd_3}
        # except:
        top_1 = {"top solution 1": 0, "top solution 2": 0, "top solution 3": 0}
    # top_1 = {"tfz": data_1["tfz"], "llg": data_1["llg"], "rmsd": rmsd_1}
    top_2 = {"top solution 1": data_1["llg"], "top solution 2": data_2["llg"], "top solution 3": data_3["llg"]}
    top_3 = {"top solution 1": data_1["tfz"], "top solution 2": data_2["tfz"], "top solution 3": data_3["tfz"]}
    # top_1 = {"top solution 1": rmsd_1, "top solution 2": rmsd_2, "top solution 3": rmsd_3}
    # top_2 = {"top solution 2": rmsd_2}
    # top_3 = {"top solution 3": rmsd_3}
    # # top_2 = {"tfz": data_2["tfz"], "llg": data_2["llg"], "rmsd": rmsd_2}
    # top_2 = {"rmsd": rmsd_2,  "tfz": data_2["tfz"]}
    # # top_3 = {"tfz": data_3["tfz"], "llg": data_3["llg"], "rmsd": rmsd_3}
    # top_3 = {"rmsd": rmsd_3,  "tfz": data_3["tfz"]}
    list_dir_1.append(top_1)
    list_dir_2.append(top_2)
    list_dir_3.append(top_3)

data_1 = pd.DataFrame(list_dir_1)
data_2 = pd.DataFrame(list_dir_2)
data_3 = pd.DataFrame(list_dir_3)

data_1.plot()
data_2.plot()
data_3.plot()
plt.show()

# print(dir_dict)
