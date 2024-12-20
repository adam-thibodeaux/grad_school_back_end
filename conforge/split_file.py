import re
from Bio.PDB import PDBIO, Structure, Model, Chain, Atom, PDBParser
import json

input_file_name = "./paritaprevir_conformers_trial_2.pdb"
output_file_path = "./para_conformers_large_trial_2/paritaprevir_conforge_"
def split_pdb_from_sdf():
    with open(input_file_name,'r') as file:
        conformers = file.read().split("END")
        file.close()

    for i in range(len(conformers)):
        with open(output_file_path + str(i) + ".pdb","w") as file:
            string = conformers[i].strip()
            #removing charge annotations from element symbols
            string = re.sub(r"\+[0-9]", "", string)
            #applying occupancy of 1 and temperature factor of 20%
            string = re.sub(r"0.00  0.00","1.00  20.00", string)
            string += "\nEND"
            file.write(string)
            file.close()
def split_native_pdb():
    char_list = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T"]
    data_list = []
    headers = ""
    with open(input_file_name,'r') as file:
        file_str = file.read()
        split_file_str = file_str.split("\n")
        headers = split_file_str[:5]
        for char in char_list:
            lines = [string for string in split_file_str if f'{char}UNK' in string]
            data_list.append(lines)
        file.close()
    for i in range(len(data_list)):
        datum = data_list[i]
        with open(output_file_path + str(i) + ".pdb", "w") as file:
            string = "\n".join(headers)
            string += "\n"
            string += "\n".join(datum)
            file.write(string)
def create_native_pdb():
    original_pdb_file = "/Users/adam/Downloads/kobe0065.pdb"
    parser = PDBParser()
    pdb_struct = parser.get_structure(id="UNK", file=original_pdb_file)
    locations = []
    with open("./temp.json", "r") as file:
        locations = json.load(file)
        file.close()
    print(len(locations))
    i = 0
    for location in locations:
        for model in pdb_struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom.set_coord(location[atom.get_name()])
        io = PDBIO()
        io.set_structure(pdb_struct)
        io.save(f"./kobe_conformers_large/kobe_conforge_{i}.pdb")
        i+= 1

    
        # for location in locations:
        #     i = 0
        #     print(location)
        #     structure = Structure.Structure("paritaprevir")
        #     model = Model.Model(0)
        #     chain = Chain.Chain("A")
        #     for atom_name, coords in location.items():
        #         atom = Atom.Atom(atom_name, coords, 0.0, 1.0, " ", atom_name, i+1, re.sub(r"[0-9]*", "", atom_name))
        #         chain.add(atom)
        #         i+=1
        #     model.add(chain)
        #     structure.add(model)
        #     io = PDBIO()
        #     io.set_structure(structure)
        #     io.save("./test.pdb")

create_native_pdb()