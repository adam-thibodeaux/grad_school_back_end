import os
import pathlib
import re
from dash import Dash, dash_table
import pandas as pd
import rmsd_calc
import ast
from rdkit.Chem import rdDistGeom, AllChem, TorsionFingerprints
import rdkit.Chem as Chem
import ast
input_struct = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb"
#input dir need trailing backslash
input_dir = "/Users/adam/Downloads/outputs_from_molec_replac/KOBE_TRIAL_1/"
output_dir = "/Users/adam/Downloads/outputs_from_molec_replac/KOBE_TRIAL_1/summary"
tfd_file_path="/Users/adam/Downloads/inputs_for_molec_replac/PAR_TFD_TRIAL_1/tfd_out.txt"

def extract_info(string, ensemble_name, iteration, solu_dict):
    string = string.strip()
    if len(string) < 5: return
    try:
        result = re.findall(r"TFZ==[0-9]*\.[0-9]*", string)
        solu_dict[ensemble_name] = solu_dict.get(ensemble_name, {"tfz": [], "llg": [], "space_group": [], })
        tfz_score = 0
        for res in result:
            tfz_score = max(float(res.split("==")[-1]), tfz_score)
        solu_dict[ensemble_name]["tfz"].append((iteration, tfz_score))
        result = re.findall(r"LLG=[0-9]*\.?[0-9]*", string)
        llg_score = -100000
        for res in result:
            llg_score = max(float(res.split("=")[-1]), llg_score)
        solu_dict[ensemble_name]["llg"].append((iteration, llg_score))
        result = re.findall(r"SOLU SPAC .*", string)

        for res in result:
            res = re.sub(r"SOLU SPAC|SOLU.*", "", res)
            solu_dict[ensemble_name]["space_group"].append((iteration, res.strip()))
        if len(solu_dict[ensemble_name]["tfz"]) > 5:
            solu_dict[ensemble_name]["tfz"] = sorted(solu_dict[ensemble_name]["tfz"], key = lambda x : x[-1], reverse=True)[:5]
            valid_iterat = []
            for tfz in solu_dict[ensemble_name]["tfz"]:
                iterat = tfz[0]
                valid_iterat.append(iterat)
            solu_dict[ensemble_name]["llg"] = [llg for llg in solu_dict[ensemble_name]["llg"] if llg[0] in valid_iterat]
            solu_dict[ensemble_name]["space_group"] = [space_group for space_group in solu_dict[ensemble_name]["space_group"] if space_group[0] in valid_iterat]

        return
    except Exception as e:
        return 

    
def main(input_struct=input_struct, input_dir=input_dir, output_dir=output_dir, should_calc_rmsd=False, is_serving=False, should_find_tfd=False, tfd_file_path=tfd_file_path):
    tfd_vals = {}
    with open(tfd_file_path) as file:
        temps = ast.literal_eval(file.read())
        for temp in temps:
            tfd_vals[temp[0]] = temp[1]
        file.close()
    solu_dict = {}
    temp_solu_dict = {}
    print(input_dir)
    for f in os.listdir(input_dir):
        if f[0] == "." or "summary" in f: continue
        for inner_f in os.listdir(input_dir + f):
            suff = pathlib.Path(inner_f).suffix
            if suff != ".sol": continue
            with open(input_dir + f + "/" + inner_f) as file:
                sol = file.read()
                sol = re.sub(r"# \[.*\]|SOLU RESOLUTION [0-9].[0-9].*|CLUSTER [0-9].*|\n", "", sol)
                sol = re.sub(r"#TFZ", "TFZ", sol)
                sol = sol.strip()
                solu_list = sol.split("#")
                temp_solu_dict[inner_f] = solu_list
        

    # sorted_solu = sorted(solu_dict.values(), key = lambda x : find_tfz(x[0]))
    for filename, solu in temp_solu_dict.items():
        for i in range(0,5):
            try:
                inner_solu = solu[i]
                result = extract_info(inner_solu, filename.replace(".sol", ""), i+ 1, solu_dict=solu_dict)
            except:
                continue

    flat_solu_dict = []
    for solu_name, solu_val in solu_dict.items():    
        for tfz in solu_val["tfz"]:
            temp_flat = {}
            iterat = tfz[0]
            temp_flat["name"] = solu_name + "." + str(iterat)
            temp_flat["tfz"] = tfz[1]
            if should_calc_rmsd:
                temp_flat["rmsd"] = rmsd_calc.calc_rmsd(input_struct_dir=input_struct, compare_struct_dir=input_dir+solu_name.split("_out")[0] + "/" + solu_name+"."+str(iterat + 1) + ".pdb")
            if should_find_tfd:
                temp_flat["tfd"] = tfd_vals.get(int(solu_name.split("_")[-2]),1)
                # print(solu_name.split("_"))[-2]
            try:
                match_llg = [llg for llg in solu_val["llg"] if llg[0] == iterat][0]
                temp_flat["llg"] = match_llg[1]
            except:
                ""
            try:
                match_space_group = [space_group for space_group in solu_val["space_group"] if space_group[0] == iterat][0]
                temp_flat["space_group"] = match_space_group[1]
            except:
                ""
            flat_solu_dict.append(temp_flat)
    df = pd.DataFrame(flat_solu_dict)
    df.sort_values("llg", inplace=True, ascending=False)
    if is_serving:
        app = Dash(__name__)
        app.layout = dash_table.DataTable(df.to_dict('records'), style_cell={'textAlign': 'left',
        'width': '{}%'.format(len(df.columns)),
        'textOverflow': 'ellipsis',
        'overflow': 'hidden'})
        app.run(debug=True)
    return df
# print(df.head())

def main_nested(input_struct=input_struct, input_dir=input_dir, output_dir=output_dir, should_calc_rmsd=False, is_serving=False):
    solu_dict = {}
    temp_solu_dict = {}
    print(input_dir)
    for upper_f in os.listdir(input_dir):
        if upper_f[0] == ".": continue
        for f in os.listdir(input_dir + upper_f):
            if f[0] == ".": continue
            for inner_f in os.listdir(input_dir + upper_f + "/"+  f):
                suff = pathlib.Path(inner_f).suffix
                if suff != ".sol": continue
                try:
                    with open(input_dir + upper_f + "/" + f +"/" + inner_f) as file:
                        sol = file.read()
                        sol = re.sub(r"# \[.*\]|SOLU RESOLUTION [0-9].[0-9].*|CLUSTER [0-9].*|\n", "", sol)
                        sol = re.sub(r"#TFZ", "TFZ", sol)
                        sol = sol.strip()
                        solu_list = sol.split("#")
                        temp_solu_dict[inner_f] = solu_list
                except:
                    ""  
            

        # sorted_solu = sorted(solu_dict.values(), key = lambda x : find_tfz(x[0]))
        for filename, solu in temp_solu_dict.items():
            for i in range(5):
                try:
                    inner_solu = solu[i]
                    result = extract_info(inner_solu, filename.replace(".sol", ""), i, solu_dict=solu_dict)
                except:
                    continue

        flat_solu_dict = []
        for solu_name, solu_val in solu_dict.items():    
            for tfz in solu_val["tfz"]:
                temp_flat = {}
                iterat = tfz[0]
                temp_flat["name"] = solu_name + "." + str(iterat)
                temp_flat["tfz"] = tfz[1]
                if should_calc_rmsd:
                    # temp_flat["rmsd"] = rmsd_calc.calc_rmsd(input_struct_dir=input_struct, compare_struct_dir=input_dir+solu_name.split("_out")[0] + "/" + solu_name+"."+str(iterat + 1) + ".pdb")
                    ""
                try:
                    match_llg = [llg for llg in solu_val["llg"] if llg[0] == iterat][0]
                    temp_flat["llg"] = match_llg[1]
                except:
                    ""
                try:
                    match_space_group = [space_group for space_group in solu_val["space_group"] if space_group[0] == iterat][0]
                    temp_flat["space_group"] = match_space_group[1]
                except:
                    ""
                flat_solu_dict.append(temp_flat)
    df = pd.DataFrame(flat_solu_dict)
    df.sort_values("llg", inplace=True, ascending=False)
    if is_serving:
        app = Dash(__name__)
        app.layout = dash_table.DataTable(df.to_dict('records'), style_cell={'textAlign': 'left',
        'width': '{}%'.format(len(df.columns)),
        'textOverflow': 'ellipsis',
        'overflow': 'hidden'})
        app.run(debug=True)
    return df
# print(df.head())

def extract_rmsd_only(compare_struct_path, input_dir, should_calc_tfd=False, should_reset=False, is_serving=False, sort_label="rmsd"):

    print('starting extraction')
    print(compare_struct_path)
    # return
    data = []
    # mol_template = None
    # compare_struct = None
    # if should_calc_tfd:
        # compare_struct = Chem.MolFromPDBFile(compare_struct_path, removeHs=False, proximityBonding=False)
        # mol_template = Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds.mol", sanitize=False, removeHs=False)
        # compare_struct = AllChem.AssignBondOrdersFromTemplate(mol_template, compare_struct)
    try:
        if should_calc_tfd and not should_reset:
            with open(input_dir + "data_rmsd_tfd.txt", "r") as file:
                data = ast.literal_eval(file.read())
                file.close()
        elif not should_reset:
            with open(input_dir + "data_rmsd.txt", "r") as file:
                data = ast.literal_eval(file.read())
                file.close()
    except:
        pass
    if len(data) == 0:
        for f in os.listdir(input_dir):
            file_path = input_dir + f
            try:
                if should_calc_tfd:
                    input_struct = Chem.MolFromPDBFile(file_path, removeHs=False, proximityBonding=False, sanitize=False)
                    # mol_template = Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds.mol", sanitize=False, removeHs=False)
                    compare_struct = Chem.MolFromPDBFile(compare_struct_path, removeHs=False, proximityBonding=False, sanitize=False)
                    # compare_struct = AllChem.AssignBondOrdersFromTemplate(mol_template, compare_struct)
                    # input_struct = AllChem.AssignBondOrdersFromTemplate(mol_template, input_struct)
                    tfd = TorsionFingerprints.GetBestTFDBetweenMolecules(compare_struct,input_struct, useWeights=False)
                    data.append({"name": f, "rmsd": rmsd_calc.calc_rmsd(compare_struct_path, file_path), "tfd": tfd})
                else:
                    data.append({"name": f, "rmsd": rmsd_calc.calc_rmsd(compare_struct_path, file_path)})
            except Exception as e:
                print(e)
                print(f'failed for {f}')
    df = pd.DataFrame(data)
    # print(df.head())
    df.sort_values(sort_label, ascending=True, inplace=True)
    
    if should_calc_tfd:
        with open(input_dir + "data_rmsd_tfd.txt", "w") as file:
            file.write(str(data))
            file.close()
    else:
        with open(input_dir + "data_rmsd.txt", "w") as file:
            file.write(str(data))
            file.close()
    if not is_serving:
        return df
    app = Dash(__name__)
    app.layout = dash_table.DataTable(df.to_dict('records'), style_cell={'textAlign': 'left',
        'width': '{}%'.format(len(df.columns)),
        'textOverflow': 'ellipsis',
        'overflow': 'hidden'})
    app.run(debug=True)

if __name__ == "__main__":
    main(should_calc_rmsd=False, is_serving=True, input_dir="/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CONFORGE/", should_find_tfd=False)
    # extract_rmsd_only("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb", "/Users/adam/Downloads/inputs_for_molec_replac/PAR_CUSTOM_CONF_TRIAL_5/ROUND_10/", should_calc_tfd=True, should_reset=False, is_serving=True, sort_label="tfd")
    
