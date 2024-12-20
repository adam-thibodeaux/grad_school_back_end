import gemmi
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem, rdForceFieldHelpers, rdchem, rdmolops
from rdkit.ForceField import rdForceField
import extract_sols as extract
import matplotlib.pyplot as plt
import numpy as np
import time
import math
import remove_model
import extract_sols
import os
import remove_model
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

# def calculate_intensity(map,x,y,z):
    
allowable_dihedral_angles = [
   
     #ninth one, pyrole ring
    # ["O1","C17","N1","C1"],
    
     # tenth one, planarization of pyridine ring
    ["C1","C2","O7","C28"],
    #eighth one, pyridine ring
    ["N1","C1","C2","O7" ],
     #first one, this is for the pyrimdine side chain as a while
    ["C15","C16", "N4", "C18"],
    #second one,  this is for the pyrimidine side chain amide bond angle
    ["O6","C18", "C19", "C20"],

     #fourth one, this is for the sulfone side chain
    ["N2","C6","C24","N3"],
    #fourth one part 2, this is the cyclopropane ring on sulfone side chain
    ["C24","N3","S1","C25"],
    #fifth one, this is for the amide link in macrocycle
    ["C8","C6","N2","C5"],
     
    # ["C16", "N4", "C18", "C19"],
    #third one, this is for the planarization of pyrimidine side chain
    # ["C15","C16", "N4", "C18"],
    #    ["N5","C19","C18","N4"],
   
    # #fourth one part 3, this is redoing sulfone side chain
    # ["N2","C6","C24","N3"],
    # ["N1","C17","C16", "N4"],
   
    
     #seventh one, second part of distal side of macrocycle core
     ["C17","C16","C15","C14"],
    # ["C11","C12","C13","C14"],
    # sixth one, distal side of macrocycle
    # ["C8","C9","C10","C11"],
    
    ["C15","C14","C13", "C12"],   
    ["C10","C11","C12","C13"],
     #ninth one, pyrole ring, again
    ["O1","C17","N1","C1"],
      #seventh one, second part of distal side of macrocycle core
    # ["C10","C11","C12","C13"],
   
    # ["O2","C5","N2","C6"]
    # ["C3","C2","O7", "C28"],
    # ["C1", "C2", "O7", "C28"]
    # ["C13","C14","C15","C16"]
    # ["C5","N2","C6","C24"]
    # ["C8","C9","C10","C11"]
    # ["H41", "N2","C5","O2"]
    # ["N1", "C4","C5", "O2" ]
    # ["H41", "N2","C6","C24"]
    # ["N1", "C4", "C5", "N2"]
    # ["C4", "N1", "C17", "C16"]
    # ["C13", "C14", "C15", "C16"]
    # ["C4","N1","C1","C2", ]
    # ["N1","C1","C2","O7" ]
    # ["C3","C2", "C1", "N1"]
# ["C5","N2","C6","C7"]
# ["N1","C1","C2","C17" ]
# ["C9","C10","C11","C12"]

                             ]
split_mol_indices = [
    # tenth one, planarization of pyridine ring
    [],
    #eighth one, pyridine ring
    [["C3","C4"]],
    #first one, this is for the pyrimidine side chain amide bond angle
    [["C15", "C14"]],
    ##second one, this is for the pyrimdine side chain as a whole
    [],
    #fourth one, this is for the sulfone side chain
    [["N2","C5"]],
    #fourth one part 2, this is the cyclopropane ring on sulfone side chain
    [],
    #fifth one, this is for the amide link in macrocycle core
    [["C5","C4"]],
    #ninth one, pyrole ring
    # [["O7","C2"],["C4","C5"]],
    #third one, this is for the planarization of pyrimidine side chain
    # [["C15","C14"]],
    # [],
    
    
     
    #fourth one part 3, more planarization of pyrimidine side chain
    # [["C15", "C14"]],
    
   
    #seventh one, second distal part of macrocycle core
    [["C14","C13"]],
    #sixth one, distal macrocycle core
    [["C12","C11"]],
    [["C13","C14"]],
    #ninth one, pyrole ring, again
    [["O7","C2"],["C4","C5"]],
     
    #seventh one, second distal part of macrocycle core
    # [["C15","C14"]],
   
    
    
    
    ]
set_dihedral_angles = [
    ["C13", "C14", "C15","C16"],
    ["C5","N2","C6","C24"]
                       ]
allowed_neighborhoods = [
    ["O7","C28","C29","C30","C31","C32","C33","C34","C35","C36","C37","C38","C39","C40", "N7"],
    ["N4", "O6", "N6", "C18","C19","C20","C21","C22","C23"],
    ["C24","C25","C26","S1","O4","O5", "N3"]
]
all_conformers = []
input_mtz_path = "/Users/adam/Downloads/outputs_from_molec_replac/PAR_CONFORGE_TRIAL_4/paritaprevir_conforge_4/paritaprevir_conforge_4_out.1.mtz"
# input_struct_path = "/Users/adam/Downloads/outputs_from_molec_replac/PAR_CONFORGE_TRIAL_4/paritaprevir_conforge_244/paritaprevir_conforge_244_out.1.pdb"
reference_struct_path = "/Users/adam/Downloads/paritaprevir_correct_bonds.mol"
#needs trailing backslash
output_conf_dir = "/Users/adam/Downloads/inputs_for_molec_replac/PAR_CUSTOM_CONF_TRIAL_5/ROUND_1/"
output_phaser_dir = "/Users/adam/Downloads/outputs_from_molec_replac/PAR_CUSTOM_CONF_TRIAL_4/"
# temp = []
# for neighborhood in allowed_neighborhoods:
#     for neighbor in neighborhood:
#         if "-" in neighborhood:
#             print('expanding ' + neighbor )
            
# new_conformers =[]
# for angle in allowable_dihedral_angles:
#     new_set = set()
#     for name in angle:
#         new_set.add(name)
#     temp.append(new_set)
#     new_set = set()
# allowable_dihedral_angles = temp

def setDihedralForOneBranch(rw_mol, conf, atom_indices, angle):
    """ 
    Set the dihedral angle for a single branch in a molecule.
    """
    removed_bonds = []
    for bnd in rw_mol.GetAtomWithIdx(atom_indices[-2]).GetBonds():
        aids = bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx()
        if aids[0] in atom_indices and aids[1] in atom_indices:
            continue
        rw_mol.RemoveBond(*aids)
        removed_bonds.append(aids + (bnd.GetBondType(),))
    # don't need the SSSR, just to know if bonds are in rings:
    Chem.FastFindRings(rw_mol)
    rdMolTransforms.SetDihedralDeg(conf, *atom_indices, angle)
    for bnd in removed_bonds:
        rw_mol.AddBond(*bnd)
def main(dihedral_angle_atom_names, split_idx_atom_names, count, perturb_angle, input_struct_path, should_proximity_bond=False, should_sanitize=True):
    # base_mol = Chem.MolFromPDBFile("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb")
    # base_mol = Chem.RemoveHs(base_mol)
    # bonds = base_mol.GetBonds()
    # print(len(bonds))
    # bond_indices = []
    # for bond in bonds:
    #     bond_indices.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    # dihedral_angles = set()
    # for i in range(len(bond_indices)):
    #     (start_1, end_1) = bond_indices[i]
    #     for j in range(i, len(bond_indices)):
    #         (start_2, end_2) = bond_indices[j]
    #         dihedral_angle = None
    #         if end_1 == start_2:
    #             dihedral_angle = (start_1,end_1,end_2)
    #             #s1-e1-e2
    #         elif start_1 == end_2:
    #             dihedral_angle =(end_1, start_1, start_2)
    #             #e1-s1-s2
    #         elif start_1 == start_2:
    #             dihedral_angle =(end_1, start_1, end_2)
    #             #e1-s1-e2
    #         elif end_1 == end_2:
    #             dihedral_angle =(start_1, end_2, start_2)
    #             #s1-e2-s2
    #         if dihedral_angle:
    #             if dihedral_angle[0] == dihedral_angle[1] or dihedral_angle[0] == dihedral_angle[2] or dihedral_angle[1] == dihedral_angle[2]:
    #                 continue
    #             dihedral_angles.add(dihedral_angle)

    # print(dihedral_angles)
    # print(len(dihedral_angles))

    # for dihedral_angle in dihedral_angles:
    #     (ang_1, ang_2, ang_3) = dihedral_angle
    #     print(rdMolTransforms.GetAngleDeg(base_mol.GetConformer(), ang_1, ang_2, ang_3))
        

    # def perturb_neighbors(angle, mol):
    

    # mtz = gemmi.read_mtz_file(input_mtz_path)
    # size = mtz.get_size_for_hkl()
    # # grid = mtz.get_f_phi_on_grid('DELFWT', "PHDELWT", size)
    # # complex_map = np.fft.ifftn(grid.array.conj())
    # # scale_factor = complex_map.size / grid.unit_cell.volume
    # # real_map = np.real(complex_map) * scale_factor
    # ccp4 = gemmi.Ccp4Map()
    # ccp4.grid = mtz.transform_f_phi_to_map("DELFWT", "PHDELWT", size)
    # # ccp4.update_ccp4_header()
    # # ccp4.write_ccp4_map("/Users/adam/Downloads/test.ccp4")
    # # print(ccp4.grid.get_nearest_point(gemmi.Position(10.0,10.0,10.0)))
    # # fig, ax = pyplot.subplots()
    # # x = np.linspace(0, ccp4.grid.unit_cell.a * 2, num=arr.shape[0], endpoint=False)
    # # y = np.linspace(0, ccp4.grid.unit_cell.b * 2, num=arr.shape[1], endpoint=False)
    # # X, Y = np.meshgrid(x, y, indexing='ij')
    # # cont = pyplot.contour(X, Y, arr[:,:,40])
    # # fig.gca().set_aspect('equal', adjustable='box')
    # # ax.clabel(cont, fontsize=10)
    # # pyplot.show()
    # print(input_struct_path)
    base_mol = Chem.MolFromPDBFile(input_struct_path, removeHs=False, sanitize=should_sanitize, proximityBonding=should_proximity_bond)
    # reference_mol = Chem.MolFromMolFile(reference_struct_path, removeHs=False)
    # base_mol = AllChem.AssignBondOrdersFromTemplate(reference_mol, base_mol)
    # print(base_mol)
    # # base_mol = Chem.RemoveHs(base_mol)
    # intensities = []
    # dihedral_angles = []
    atom_idx_dict = {}
    for i, atom in enumerate(base_mol.GetAtoms()):
        atom_name = atom.GetMonomerInfo().GetName().strip()
        atom_idx_dict[atom_name] = atom.GetIdx()
    # # for i, atom in enumerate(base_mol.GetAtoms()):
    # #     # if atom.GetSymbol() == "H": continue
        
    # #     # print(atom.GetMonomerInfo().GetName())
    # #     # print(atom.GetIdx())
    # #     atom_name = atom.GetMonomerInfo().GetName().strip()
    # #     atom_idx_dict[atom_name] = atom.GetIdx()
    # #     dihedral_angle = set()
    # #     # dihedral_angle.add(atom_name)
    # #     for first_neighbor in atom.GetNeighbors():
    # #         # if atom.GetSymbol() == "H": continue
    # #         # print(first_neighbor.GetMonomerInfo().GetName())
    # #         first_neighbor_name = first_neighbor.GetMonomerInfo().GetName().strip()
    # #         # dihedral_angle.add(first_neighbor_name)
    # #         for second_neighbor in first_neighbor.GetNeighbors():
    # #             # if atom.GetSymbol() == "H": continue
    # #             second_neighbor_name = second_neighbor.GetMonomerInfo().GetName().strip()
    # #             # dihedral_angle.add(second_neighbor_name)
    # #             for third_neighbor in second_neighbor.GetNeighbors():
    # #                 # if atom.GetSymbol() == "H": continue
    # #                 dihedral_angle = []
    # #                 dihedral_angle.append(atom_name)
    # #                 dihedral_angle.append(first_neighbor_name)
    # #                 dihedral_angle.append(second_neighbor_name)
                    
    # #                 third_neighbor_name = third_neighbor.GetMonomerInfo().GetName().strip()
    # #                 dihedral_angle.append(third_neighbor_name)
    # #                 if len(dihedral_angle) == 4:
    # #                     dihedral_angles.append(dihedral_angle)

                


    # #     # positions = base_mol.GetConformer().GetAtomPosition(i)
    # #     # position = gemmi.Position(positions.x, positions.y, positions.z)
    
    # #         # intensity = ccp4.grid.get_value(math.floor(positions.x), math.floor(positions.y), math.floor(positions.z))
    # #     # intensity = ccp4.grid.interpolate_value(position)
    # #         # neighborhood = []
    # #         # neighborhood.append([math.floor(positions.x), math.floor(positions.y), math.floor(positions.z)])
    # #         # neighborhood.append([math.floor(positions.x), math.ceil(positions.y), math.floor(positions.z)])
    # #         # neighborhood.append([math.floor(positions.x), math.floor(positions.y), math.ceil(positions.z)])
    # #         # neighborhood.append([math.ceil(positions.x), math.floor(positions.y), math.floor(positions.z)])
    # #         # neighborhood.append([math.ceil(positions.x), math.ceil(positions.y), math.floor(positions.z)])
    # #         # neighborhood.append([math.ceil(positions.x), math.floor(positions.y), math.ceil(positions.z)])
    # #         # neighborhood.append([math.floor(positions.x), math.ceil(positions.y), math.ceil(positions.z)])
    # #         # neighborhood.append([math.ceil(positions.x), math.ceil(positions.y), math.ceil(positions.z)])
    # #         # neighborhood_intensities = []
    # #         # for point in neighborhood:
    # #         #     neighborhood_intensities.append(ccp4.grid.get_value(*point))
    # #         #     break
    # #         # print(neighborhood_intensities)
    # #         # avg_intensity = sum(neighborhood_intensities)/len(neighborhood_intensities)
    # #         # print(avg_intensity)

    # #     # intensities.append((f"{atom.GetSymbol()} at [{positions.x, positions.y, positions.z}]", intensity))
            
    # # for el in sorted(intensities, key=lambda x: x[1]):
    #     # print(el)
    # found_dihedrals = []
    # found_set_dihedrals = []
    # for i in range(len(allowable_dihedral_angles)):
    #     dihedral = allowable_dihedral_angles[i]
    #     try:
    #         # print(f"matches for {dihedral}")
    #         matching_dihedrals = [found_dihedral for found_dihedral in dihedral_angles if len(intersection(dihedral,found_dihedral)) == 4]
    #         # print(matching_dihedrals)
    #         matching_dihedral = []
    #         print(matching_dihedrals)
    #         for name in matching_dihedrals[0]:
    #             # print(name)
    #             matching_dihedral.append(atom_idx_dict.get(name,-1))
    #         found_dihedrals.append(matching_dihedral)
    #     except:
    #         continue
    # for i in range(len(set_dihedral_angles)):
    #     dihedral = set_dihedral_angles[i]
    #     try:
    #         # print(f"matches for {dihedral}")
    #         matching_dihedrals = [found_dihedral for found_dihedral in dihedral_angles if len(intersection(dihedral,found_dihedral)) == 4]
    #         # print(matching_dihedrals)
    #         matching_dihedral = []
    #         print(matching_dihedrals)
    #         for name in matching_dihedrals[0]:
    #             # print(name)
    #             matching_dihedral.append(atom_idx_dict.get(name,-1))
    #         found_set_dihedrals.append(matching_dihedral)
    #     except:
    #         continue
    # # print(found_dihedrals)
    # # print(atom_idx_dict)
    dihedral_angle_atom_indices = [[atom_idx_dict[atom_name] for atom_name in dihedral_angle_atom_names]]

    # for found_dihedral in found_dihedrals:
    #     orig_dihedral_angles.append(rdMolTransforms.GetDihedralDeg(base_mol.GetConformer(), *found_dihedral))
    for i in range(len(dihedral_angle_atom_indices)):
        # angle = orig_dihedral_angles[i]
        angle_indices = dihedral_angle_atom_indices[i]
        # print(angle)
        # atom_indices = [atom_idx_dict[x] for x in allowable_dihedral_angles[0]]
        comp_id_pos = base_mol.AddConformer(base_mol.GetConformer(0), assignId=True)
        # comp_id_neg = base_mol.AddConformer(base_mol.GetConformer(0), assignId=True)
        edit_mol = rdchem.RWMol(base_mol)
       
        
        try:
            # setDihedralForOneBranch(edit_mol, edit_mol.GetConformer(comp_id_pos), atom_indices, perturb_angle)
            rdMolTransforms.SetDihedralDeg(base_mol.GetConformer(comp_id_pos),*angle_indices, perturb_angle)
            Chem.MolToPDBFile(base_mol, f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count}.pdb", confId=comp_id_pos)
            return base_mol
        except Exception as e:
            print("failed setting explicit dihedral angle, using MMFF optimization")
            # molec_props = rdForceFieldHelpers.MMFFGetMoleculeProperties(base_mol)
            # force_field = rdForceFieldHelpers.MMFFGetMoleculeForceField(base_mol, molec_props, confId=comp_id_pos)
            # force_field.MMFFAddTorsionConstraint(*angle_indices, relative=True, minDihedralDeg= perturb_angle, maxDihedralDeg= perturb_angle, forceConstant=0.001)
            # for j in range(len(found_set_dihedrals)):
            #     set_angle_indices = found_set_dihedrals[j]
            #     print(j, set_angle_indices)
            #     force_field.MMFFAddTorsionConstraint(*set_angle_indices, relative=True, minDihedralDeg=0, maxDihedralDeg=0, forceConstant=0.001)
            
            for split_mol_idx in split_idx_atom_names:
                bi = atom_idx_dict[split_mol_idx[0]]
                ei = atom_idx_dict[split_mol_idx[1]]

                edit_mol.RemoveBond(bi,ei)
            # return edit_mol
                ba = edit_mol.GetAtomWithIdx(bi)
                ea = edit_mol.GetAtomWithIdx(ei)
            
                ba.SetNumExplicitHs(ba.GetNumExplicitHs() + 1)
                ea.SetNumExplicitHs(ea.GetNumExplicitHs() + 1)
                Chem.SanitizeMol(edit_mol)
            setDihedralForOneBranch(edit_mol, edit_mol.GetConformer(comp_id_pos), angle_indices, perturb_angle)
            # rdMolTransforms.SetDihedralDeg(edit_mol.GetConformer(comp_id_pos),*angle_indices,perturb_angle)
            for split_mol_idx in split_idx_atom_names:

                bi = atom_idx_dict[split_mol_idx[0]]
                ei = atom_idx_dict[split_mol_idx[1]]
            # return edit_mol
                ba = edit_mol.GetAtomWithIdx(bi)
                ea = edit_mol.GetAtomWithIdx(ei)
            
                ba.SetNumExplicitHs(ba.GetNumExplicitHs() - 1)
                ea.SetNumExplicitHs(ea.GetNumExplicitHs() - 1)
                edit_mol.AddBond(bi, ei, Chem.BondType.SINGLE)
                conf = edit_mol.GetConformer(comp_id_pos)
                for bond in edit_mol.GetBonds():
                    start = bond.GetBeginAtomIdx()
                    end = bond.GetEndAtomIdx()
                    start_coords = np.array(conf.GetAtomPosition(start))
                    end_coords = np.array(conf.GetAtomPosition(end))
                    dist = np.linalg.norm(end_coords-start_coords)
                    if dist > 3:
                        print(f'suspicious bond distance {dist} angstroms for: paritaprevir_torsion_angle_perturb_{count}.pdb')
                        # print('skipping this one')
                        # return
                # Chem.SanitizeMol(edit_mol, )
            # edit_mol = AllChem.AssignBondOrdersFromTemplate(reference_mol, edit_mol)
            

            # edit_mol.GetConformer(comp_id_pos).SetAtomPosition(bi,ba_coords)
            # edit_mol.GetConformer(comp_id_pos).SetAtomPosition(ei,ea_coords)
            # print(ba_coords.x)
            # print(ea_coords.x)
            Chem.SanitizeMol(edit_mol)
            rdForceFieldHelpers.MMFFOptimizeMolecule(edit_mol)
            
            Chem.MolToPDBFile(edit_mol, f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count}.pdb", confId=comp_id_pos)
            return edit_mol


        # rdMolTransforms.SetDihedralDeg(base_mol.GetConformer(comp_id_neg),*angle_indices,angle - perturb_angle)
    # count = 0
    # total_count = 0
    # for i in range(len(dihedral_angle_atom_indices)):
    #     # print(f"original angle {orig_dihedral_angles[i]}")
    #     angle_indices = found_dihedrals[i]
    #     # for j in range(base_mol.GetNumConformers()):
    #         # print(rdMolTransforms.GetDihedralDeg(base_mol.GetConformer(j),*angle_indices))
    # total_base_conformer_intensity = []
    # total_new_conformer_intensity = []
    # return base_mol
    # for j in range(base_mol.GetNumConformers()):
    #     total_base_conformer_intensity.append([])
    #     total_new_conformer_intensity.append([])
    #     for neighborhood in allowed_neighborhoods:
    #         base_neighborhood_intensities = []
    #         for neighbor in neighborhood:
    #             try:
    #                 atom_idx = atom_idx_dict[neighbor]
    #                 base_position = base_mol.GetConformer(0).GetAtomPosition(atom_idx)
    #                 base_intensity = ccp4.grid.interpolate_value(gemmi.Position(base_position.x,base_position.y,base_position.z))
    #                 total_base_conformer_intensity[j].append(base_intensity)
    #                 base_neighborhood_intensities.append(base_intensity)
                
    #                 neighborhood_intensities = []
    #                 position = base_mol.GetConformer(j).GetAtomPosition(atom_idx)
    #                 # print(position.x - base_position.x, position.y-base_position.y, position.z-base_position.z)
    #                 intensity = ccp4.grid.interpolate_value(gemmi.Position(position.x,position.y,position.z))
    #                 total_count += 1
    #                 neighborhood_intensities.append(intensity)
    #                 total_new_conformer_intensity[j].append(intensity)
    #                 if intensity > base_intensity:
    #                     # print("i think this is better")
    #                     # print(f"{neighbor} in conformer {j}")
    #                     # print(f"previous intensity of {base_intensity} to {intensity}")
    #                     count += 1
        
    #             except:
    #                 print("errored")
    #                 continue
    #         avg_base_intensity = sum(base_neighborhood_intensities)/len(base_neighborhood_intensities)
    #         avg_new_intensity = sum(neighborhood_intensities)/len(neighborhood_intensities)
    #         # print(f"neighborhood intensity for {neighborhood} at conformer {j}, changed from {avg_base_intensity} to {avg_new_intensity} ")
    #         # print(f"this is a {((avg_new_intensity - avg_base_intensity)/avg_base_intensity)*100}% change")
    #         # if avg_new_intensity > avg_base_intensity:
    #             # print("good change")
    #         # else:
    #             # print("bad change")
    # print(f"{count} of {total_count} atoms improved intensity") 
    # total_conformer_intensity_change = []
    # print(f"summary for all conformers for perturb angle of {perturb_angle} degrees")
    # for i in range(len(total_base_conformer_intensity)):
    #     base_conformer_intensity = total_base_conformer_intensity[i]
    #     new_conformer_intensity = total_new_conformer_intensity[i]
    #     base_intensity_sum = sum(base_conformer_intensity)
    #     new_intensity_sum = sum(new_conformer_intensity)
    #     print(f"conformer #{i} changed from {base_intensity_sum} to {new_intensity_sum} intensity")
    #     total_conformer_intensity_change.append(new_intensity_sum-base_intensity_sum)
    # all_conformers.append(base_mol)
    # return total_conformer_intensity_change
    #         # writer = Chem.PDBWriter(f"/Users/adam/Downloads/inputs_for_molec_replac/custom_confs/paritaprevir_custom_{j}.pdb")
    #         # writer.write(mol=base_mol, confId=j)
        
    
    ##convert this to np array
    ##check for local minima/maxima
    ##do coords match pdb file?
    ## perturb -> how much to perturb?
    ##rerun through phaser

def generate_conformers(count, step, dihedral_angle_atom_names, split_idx_atom_names, input_struct_path, should_proximity_bond=False):
    for run_num in range(0,360,step):
        count += 1
        try:

            main(dihedral_angle_atom_names=dihedral_angle_atom_names, split_idx_atom_names=split_idx_atom_names, count=count, perturb_angle=run_num, input_struct_path=input_struct_path, should_proximity_bond=should_proximity_bond)

            # Chem.MolToPDBFile(mol, f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count}.pdb", confId=1)
            continue
        except Exception as e:
            print(e)
            continue
    return count

def extract_rmsd_and_update_init_path(input_struct_path):
    input_struct_path_list =set()
    data = extract.extract_rmsd_only(input_struct_path,output_conf_dir, should_calc_tfd=False, should_reset=True)
    print(data.to_string())
    var = np.var(data["rmsd"])
    vals =[np.mean(data["rmsd"])]
    for i in range(1,6):
        vals.append(vals[0] - i * var)
        vals.append(vals[0] + i * var)
    for val in vals:
        input_struct_path_list.add(output_conf_dir +data.iloc[(data['rmsd']-val).abs().argsort()[:2]]["name"].tolist()[0])
    print(input_struct_path_list)
    
    plt.hist(data["rmsd"], bins=36)
    plt.show()
    return input_struct_path_list
def extract_llg_and_tfz_and_update_init_path(input_dir, num_to_include=1):
    data = extract.main(input_dir=input_dir)
    data["tfz_cross_llg"] = data["tfz"] * data["llg"]/100
    
    with open(input_dir + "summary.txt", "w") as f:
        f.write(data.to_string())
        f.close()
    # data = data[data["tfz"] >= 5]
    data.sort_values("llg", ascending=False, inplace=True)
    print(data.to_string())
    set_input_path_list = set()
    i = 0
    while len(set_input_path_list) < num_to_include and i < len(data):
        set_input_path_list.add(data.iloc[i]["name"].split("_out")[0])
        i += 1
    file_prefix = input_dir.split("PHASER")[0]
    # return [f{file_prefix}]
    return [f"{file_prefix}{x}.pdb" for x in set_input_path_list]
    ret_list = [input_dir + data.iloc[0]["name"].split("_out")[0] +"/"+ data.iloc[0]["name"]+".pdb",
            input_dir + data.iloc[1]["name"].split("_out")[0] +"/"+ data.iloc[1]["name"]+".pdb",
            input_dir + data.iloc[2]["name"].split("_out")[0] +"/"+ data.iloc[2]["name"]+".pdb"]
    return ret_list[:num_to_include]
   
input_struct_path_list = []
for i in range(len(allowable_dihedral_angles)):
    
    input_struct_path =  "/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CONFORGE/paritaprevir_conforge_13/paritaprevir_conforge_13_out.1.pdb"
    input_mtz = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta.mtz"
    count = 0
    dihedral_angle = allowable_dihedral_angles[i]   
    split_idx = split_mol_indices[i]
    output_conf_dir = f"/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CUSTOM_CONF/ROUND_{i}/"
    if not os.path.isdir(output_conf_dir):
        os.mkdir(output_conf_dir)
    if not os.path.isdir(output_conf_dir + "PHASER"):
        os.mkdir(output_conf_dir + "PHASER")
    if i > 0:
        #find top5 from previous round based on tfd/rmsd, loop through with step=5
        for input_struct_path_individ in input_struct_path_list:

            input_struct_path = input_struct_path_individ
            if i <= 1:
                count = generate_conformers(count=count, step=20, dihedral_angle_atom_names=dihedral_angle,split_idx_atom_names=split_idx, input_struct_path=input_struct_path, should_proximity_bond=True)
            else:
                count = generate_conformers(count=count, step=20, dihedral_angle_atom_names=dihedral_angle,split_idx_atom_names=split_idx, input_struct_path=input_struct_path, should_proximity_bond=False)
        # input_struct_path_list = list(extract_rmsd_and_update_init_path(input_struct_path))
        input_struct_path_list = output_conf_dir
        # for j in range(len(input_struct_path_list)):
        timers = []
        sorted_output_conf = []
        for j in os.listdir(output_conf_dir):
            sorted_output_conf.append(j)
        sorted_output_conf = sorted(sorted_output_conf)
        for k in range(len(sorted_output_conf)):
            j = sorted_output_conf[k]
            phaser_input = output_conf_dir + j
            if 'paritaprevir' not in phaser_input: continue
            phaser_name = phaser_input.split("/")[-1].replace(".pdb","")
            if os.path.isdir(output_conf_dir + "PHASER/" + phaser_name):
                print(f'PHASER already done for {phaser_name}, skipping, round {i}')
                continue
            print('starting PHASER for ' + phaser_name)
            timer = remove_model.main_tfd(input_pdb=phaser_input,output_dir=output_conf_dir+"PHASER/"+phaser_name+"/", out_dir= phaser_name, input_mtz_tfd=input_mtz)
            timers.append(timer)
            print(f'finished PHASER for {phaser_name} in {timer} seconds')
            print(f'approximate time remaining for this round {i} is: {(sum(timers)/len(timers)) * (len(sorted_output_conf) - k)} seconds')
        if i <= 6:
            input_struct_path_list =extract_llg_and_tfz_and_update_init_path(output_conf_dir + "PHASER/", num_to_include=3)
        else:
            input_struct_path_list =extract_llg_and_tfz_and_update_init_path(output_conf_dir + "PHASER/", num_to_include=3)
    else:
        generate_conformers(count=count, step=20, dihedral_angle_atom_names=allowable_dihedral_angles[i], split_idx_atom_names=split_mol_indices[i], input_struct_path=input_struct_path, should_proximity_bond=True)
        # input_struct_path_list = list(extract_rmsd_and_update_init_path(input_struct_path))
        input_struct_path_list = output_conf_dir
        # for j in range(len(input_struct_path_list)):
        timers = []
        sorted_output_conf = []
        for j in os.listdir(output_conf_dir):
            sorted_output_conf.append(j)
        sorted_output_conf = sorted(sorted_output_conf)
        for k in range(len(sorted_output_conf)):
            j = sorted_output_conf[k]
            phaser_input = output_conf_dir + j
            if 'paritaprevir' not in phaser_input: continue
            phaser_name = phaser_input.split("/")[-1].replace(".pdb","")
            if os.path.isdir(output_conf_dir + "PHASER/" + phaser_name):
                print(f'PHASER already done for {phaser_name}, skipping, round {i}')
                continue
            print('starting PHASER for ' + phaser_name)
            timer = remove_model.main_tfd(input_pdb=phaser_input,output_dir=output_conf_dir+"PHASER/"+phaser_name+"/", out_dir= phaser_name)
            timers.append(timer)
            print(f'finished PHASER for {phaser_name} in {timer} seconds')
            print(f'approximate time remaining for this round {i} is: {(sum(timers)/len(timers)) * (len(sorted_output_conf) - k)} seconds')
        input_struct_path_list =extract_llg_and_tfz_and_update_init_path(output_conf_dir + "PHASER/", num_to_include=3)

        


#         if os.path.isdir(f"{output_conf_dir}RUN_{run_num}") and os.path.isdir(f"{output_phaser_dir}RUN_{run_num}"):
#             print(f"already completed {run_num}, skipping...")
#             continue
#         os.mkdir(f"{output_conf_dir}RUN_{run_num}")
#         os.mkdir(f"{output_phaser_dir}RUN_{run_num}")
#         best_conformer_deg = [(0,0)] *(len(allowable_dihedral_angles) + 1)
#         for i in range(5,6):
#             perturb_angle = i
#             intensities = main()
#             print(intensities[1])
#             for j in range(len(intensities)):
#                 best_conformer_deg[j] = (i,intensities[j])
#         print(best_conformer_deg)
#         for i  in range(1,len(best_conformer_deg)):
#             best_conf = best_conformer_deg[i]
#             if best_conf[0] == 0: continue
#             print(best_conf[0])
#             writer = Chem.PDBWriter(f"{output_conf_dir}RUN_{run_num}/run_{run_num}_conf_{i}.pdb")
#             # all_conformers[best_conf[0] - 1] = Chem.AddHs(all_conformers[best_conf[0] - 1])
#             writer.write(all_conformers[0], confId=i)
#             # print(all_conformers[best_conf[0] - 1].GetConformer(i))
#         # for i in range(len(all_conformers)):
#             # conf_tuple = all_conformers[i]
#             # if conf_tuple[0] == 0: continue
#             # deg_index = conf_tuple[0]
#         input_phaser = f"{output_conf_dir}RUN_{run_num}/"
#         output_phaser = f"{output_phaser_dir}RUN_{run_num}/"

#         ### below is for phaser automation
#         # remove_model.main(fixed_ens_dir=input_phaser,output_dir=output_phaser)
#         # print(output_phaser)
#         # df = extract_sols.main(input_dir=output_phaser,output_dir=output_phaser + "summary")
#         # print(df.head())
#         # # print(df.iloc[0]['name'].spl)
#         # top_conf = df.iloc[0]['name']
#         # off_by_one_int = int(top_conf[-1]) + 1
#         # top_conf = top_conf[:-1]
#         # top_conf += str(off_by_one_int)
#         # top_conf_base = top_conf.split("_out")[0]
#         # input_mtz_path = f"{output_phaser}{top_conf_base}/{top_conf}.mtz"
#         # input_struct_path = f"{output_phaser}{top_conf_base}/{top_conf}.pdb"
#         # print(input_struct_path)
#         # print(input_mtz_path)
#         # extract_dir = f"{}"

        

# # # after first run initialization
# # count = 0
# # for f in os.listdir("/Users/adam/Downloads/inputs_for_molec_replac/PAR_CUSTOM_CONF_TRIAL_3/ROUND_2"):
# #     if "pdb" not in f.split(".")[-1]: continue
# #     input_struct_path = f"/Users/adam/Downloads/inputs_for_molec_replac/PAR_CUSTOM_CONF_TRIAL_3/ROUND_2/{f}"

# #     for run_num in range(0,360):
# #         count += 1
# #         try:
# #             mol = main()
# #             print(mol.GetNumConformers())
# #             perturb_angle += 1
# #             Chem.MolToPDBFile(mol, f"{output_conf_dir}paritaprevir_angle_perturb_{count}.pdb", confId=1)

# #             continue
# #         except Exception as e:
# #             print(e)
#             continue
#         if os.path.isdir(f"{output_conf_dir}RUN_{run_num}") and os.path.isdir(f"{output_phaser_dir}RUN_{run_num}"):
#             print(f"already completed {run_num}, skipping...")
#             continue
#         os.mkdir(f"{output_conf_dir}RUN_{run_num}")
#         os.mkdir(f"{output_phaser_dir}RUN_{run_num}")
#         best_conformer_deg = [(0,0)] *(len(allowable_dihedral_angles) + 1)
#         for i in range(5,6):
#             perturb_angle = i
#             intensities = main()
#             print(intensities[1])
#             for j in range(len(intensities)):
#                 best_conformer_deg[j] = (i,intensities[j])
#         print(best_conformer_deg)
#         for i  in range(1,len(best_conformer_deg)):
#             best_conf = best_conformer_deg[i]
#             if best_conf[0] == 0: continue
#             print(best_conf[0])
#             writer = Chem.PDBWriter(f"{output_conf_dir}RUN_{run_num}/run_{run_num}_conf_{i}.pdb")
#             # all_conformers[best_conf[0] - 1] = Chem.AddHs(all_conformers[best_conf[0] - 1])
#             writer.write(all_conformers[0], confId=i)
#             # print(all_conformers[best_conf[0] - 1].GetConformer(i))
#         # for i in range(len(all_conformers)):
#             # conf_tuple = all_conformers[i]
#             # if conf_tuple[0] == 0: continue
#             # deg_index = conf_tuple[0]
#         input_phaser = f"{output_conf_dir}RUN_{run_num}/"
#         output_phaser = f"{output_phaser_dir}RUN_{run_num}/"

#         ### below is for phaser automation
#         # remove_model.main(fixed_ens_dir=input_phaser,output_dir=output_phaser)
#         # print(output_phaser)
#         # df = extract_sols.main(input_dir=output_phaser,output_dir=output_phaser + "summary")
#         # print(df.head())
#         # # print(df.iloc[0]['name'].spl)
#         # top_conf = df.iloc[0]['name']
#         # off_by_one_int = int(top_conf[-1]) + 1
#         # top_conf = top_conf[:-1]
#         # top_conf += str(off_by_one_int)
#         # top_conf_base = top_conf.split("_out")[0]
#         # input_mtz_path = f"{output_phaser}{top_conf_base}/{top_conf}.mtz"
#         # input_struct_path = f"{output_phaser}{top_conf_base}/{top_conf}.pdb"
#         # print(input_struct_path)
#         # print(input_mtz_path)
#         # extract_dir = f"{}"

        

