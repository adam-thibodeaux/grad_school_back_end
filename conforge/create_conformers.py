import sys
import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen
import CDPL.Biomol as Biomol
import re
import json
import CDPL.Base as Base

def generateConformationEnsembles(mol: Chem.BasicMolecule, conf_gen: ConfGen.ConformerGenerator) -> (int, int):
     """
     Generates a conformation ensemble for the argument molecule using the provided initialized ConfGen.ConformerGenerator instance.

     Parameters:
     - mol (Chem.BasicMolecule): Molecule to generate a conformation ensemble for.
     - conf_gen (ConfGen.ConformerGenerator): Instance of the ConfGen.ConformerGenerator class.

     Returns:
     - int: Status code indicating the success of the conformation ensemble generation.
     - int: Number of generated conformers.
     """
     # prepare the molecule for conformer generation
     ConfGen.prepareForConformerGeneration(mol)

     # generate the conformer ensemble
     status = conf_gen.generate(mol)
     num_confs = conf_gen.getNumConformers()
     # if sucessful, store the generated conformer ensemble as
     # per atom 3D coordinates arrays (= the way conformers are represented in CDPKit)
     if status == ConfGen.ReturnCode.SUCCESS or status == ConfGen.ReturnCode.TOO_MUCH_SYMMETRY:

         conf_gen.setConformers(mol)
     else:
         num_confs = 0

     # return status code and the number of generated conformers
     return (status, num_confs)


def generate3dConformation(mol: Chem.Molecule, struct_gen: ConfGen.StructureGenerator) -> int:
     """
     Generates a low-energy 3D structure of the argument molecule using the provided initialized ConfGen.StructureGenerator instance.

     Parameters:
     - mol (Chem.Molecule): Molecule to generate a 3D structure for.
     - struct_gen (ConfGen.StructureGenerator): Instance of the ConfGen.StructureGenerator class.

     Returns:
     - int: Status code indicating the success of the 3D structure generation.
     """
     # prepare the molecule for 3D structure generation
     ConfGen.prepareForConformerGeneration(mol)

     # generate the 3D structure
     status = struct_gen.generate(mol)

     # if sucessful, store the generated conformer ensemble as
     # per atom 3D coordinates arrays (= the way conformers are represented in CDPKit)
     if status == ConfGen.ReturnCode.SUCCESS:
         struct_gen.setCoordinates(mol)

     # return status code
     return status

max_time = 3600 # Max. allowed molecule processing time in seconds (default: 3600 sec)

mols = []  # Example list of BasicMolecules
# reader = Biomol.FilePDBMoleculeReader("/Users/adam/Downloads/kobe0065.pdb")
reader = Chem.MoleculeReader("./kobe0065.sdf")

mol = Chem.BasicMolecule()


 # dictionary mapping status codes to human readable strings
status_to_str = { ConfGen.ReturnCode.UNINITIALIZED                  : 'uninitialized',
                   ConfGen.ReturnCode.TIMEOUT                        : 'max. processing time exceeded',
                   ConfGen.ReturnCode.ABORTED                        : 'aborted',
                   ConfGen.ReturnCode.FORCEFIELD_SETUP_FAILED        : 'force field setup failed',
                   ConfGen.ReturnCode.FORCEFIELD_MINIMIZATION_FAILED : 'force field structure refinement failed',
                   ConfGen.ReturnCode.FRAGMENT_LIBRARY_NOT_SET       : 'fragment library not available',
                   ConfGen.ReturnCode.FRAGMENT_CONF_GEN_FAILED       : 'fragment conformer generation failed',
                   ConfGen.ReturnCode.FRAGMENT_CONF_GEN_TIMEOUT      : 'fragment conformer generation timeout',
                   ConfGen.ReturnCode.FRAGMENT_ALREADY_PROCESSED     : 'fragment already processed',
                   ConfGen.ReturnCode.TORSION_DRIVING_FAILED         : 'torsion driving failed',
                   ConfGen.ReturnCode.CONF_GEN_FAILED                : 'conformer generation failed' }

 # process molecules one after the other
conf_gen = ConfGen.ConformerGenerator()

max_time = 7200 # Max. allowed molecule processing time in seconds (default: 3600 sec)
min_rmsd = 0.02 # Output conformer RMSD threshold (default: 0.5)
e_window = 20 # Output conformer energy window (default: 20.0)
max_confs = 1000 # Max. output ensemble size (default: 100)

conf_gen.settings.timeout = max_time * 1000          # apply the -t argument
conf_gen.settings.minRMSD = min_rmsd                 # apply the -r argument
conf_gen.settings.energyWindow = e_window            # apply the -e argument
conf_gen.settings.maxNumOutputConformers = max_confs # apply the -n argument
# conf_gen.settings.maxRotatableBondCount = 20
# conf_gen.settings.forceFieldTypeStochastic = 0
# conf_gen.settings.forceFieldTypeSystematic = 1
print(conf_gen.settings.getMacrocycleRotorBondCountThreshold())
print(conf_gen.settings.getForceFieldTypeStochastic())
print(conf_gen.settings.getMaxRotatableBondCount())
# print(conf_gen.settings.getMinRMSD())
 # dictionary mapping status codes to human readable strings
status_to_str = { ConfGen.ReturnCode.UNINITIALIZED                  : 'uninitialized',
                   ConfGen.ReturnCode.TIMEOUT                        : 'max. processing time exceeded',
                   ConfGen.ReturnCode.ABORTED                        : 'aborted',
                   ConfGen.ReturnCode.FORCEFIELD_SETUP_FAILED        : 'force field setup failed',
                   ConfGen.ReturnCode.FORCEFIELD_MINIMIZATION_FAILED : 'force field structure refinement failed',
                   ConfGen.ReturnCode.FRAGMENT_LIBRARY_NOT_SET       : 'fragment library not available',
                   ConfGen.ReturnCode.FRAGMENT_CONF_GEN_FAILED       : 'fragment conformer generation failed',
                   ConfGen.ReturnCode.FRAGMENT_CONF_GEN_TIMEOUT      : 'fragment conformer generation timeout',
                   ConfGen.ReturnCode.FRAGMENT_ALREADY_PROCESSED     : 'fragment already processed',
                   ConfGen.ReturnCode.TORSION_DRIVING_FAILED         : 'torsion driving failed',
                   ConfGen.ReturnCode.CONF_GEN_FAILED                : 'conformer generation failed' }

try:
    while reader.read(mol):
        locations = []
        try:
            
            print('Molecule with', mol.numAtoms, 'atoms and', mol.numBonds, 'bonds')
            # compose a simple molecule identifier
            mol_id = Chem.getName(mol).strip()
            if mol_id == '':
                mol_id = '#' + 'kobe_0065' # fallback if name is empty
            else:
                mol_id = '\'%s\' (#kobe_0065)' % (mol_id)

            try:
                #     with open("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb", "r") as file:
                #         file_str_split = file.read().split("\n")
                #         for atom_name, coords in location.items():
                #             matching_index = [idx for idx, s in enumerate(file_str_split) if atom_name in s][0]
                #             matching_nums = re.findall(coord_regex, file_str_split[matching_index])
                #             temp = matching_nums[3]
                #             b_fact = matching_nums[4]
                #             template = {file_str_split[matching_index][:31]}
                #             for coord in coords:
                #                 if coord > 0:
                #                     template += f"   {coord:.3f}"
                #                 else:
                #                     template += f"  {coord:.3f}"
                #             symbol_idx = file_str_split[matching_index].find(r" [A-Z] ")
                #             if len(temp) > 3:
                #                 template += f" {temp}  {b_fact}          {re.sub(r"[0-9]*", "", atom_name)}"
                #             # print(template)
                #             file_str_split[matching_index] = template
                #             # for i in range(len(matching_coords)):
                #             #     matching_coord = matching_coords[i]
                #             #     round_coord = str("%.3f" % coords[i])
                #             #     file_str_split[matching_index] = re.sub(matching_coord, round_coord, file_str_split[matching_index])
                #         print("\n".join(file_str_split))
                #         file.close()
                # raise
                # generate 3D structure of the read molecule
                # print()
               
                # coord_regex = re.compile(r"-?[0-9]{1,3}\.[0-9]{3,4}")
                # with open("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb", "r") as file:
                #     file_str = file.read()
                #     file_str_split = file_str.split("\n")
                #     for key, val in locations.items():
                #         matching_index = [idx for idx, s in enumerate(file_str_split) if key in s][0]
                #         matching_coords = re.findall(coord_regex, file_str_split[matching_index])

                #         for i in range(len(matching_coords)):
                #             coord = str(matching_coords[i])
                #             for j in range(len(val)):
                #                 round_val = "%.3f" % val[0][i]
                #                 file_str_split[matching_index] = re.sub(coord,str(round_val), file_str_split[matching_index])
                #             print("\n".join(file_str_split))
                        
                #     file.close()
                    
                # raise
                status, num_confs = generateConformationEnsembles(mol, conf_gen)
                print(f"generated {num_confs} conformers")
                
                # with open("./temp.txt", "w") as file:
                #     file.close()
                # check for severe error reported by status code
                atoms = mol.getAtoms()
                
                for i in range(num_confs):
                    count = {}
                    sub_location = {}
                    for atom in atoms:
                        atom_name = Chem.getSymbol(atom)
                        count[atom_name] = count.get(atom_name, 0) + 1
                        coords = Chem.getConformer3DCoordinates(atom, i)
                        sub_location[f"{atom_name}{count[atom_name]}"] = [coords.getElement(0), coords.getElement(1), coords.getElement(2)]
                    locations.append(sub_location)
                
                with open("./temp.json", "w") as file:
                    json.dump(locations, file)
                    file.close()
             # output the generated 3D structure
                # if not writer.write(mol):
                #     sys.exit('Error: writing 3D structure of molecule %s failed' % mol_id)

            except Exception as e:
                print(locations)
                sys.exit('Error: 3D structure generation or output for molecule %s failed: %s' % (mol_id, str(e)))
                
        except Exception as e:
                sys.exit('Error: processing of molecule failed: ' + str(e))
except Exception as e: # handle exception raised in case of severe read errors
    sys.exit('Error: reading molecule failed: ' + str(e))


