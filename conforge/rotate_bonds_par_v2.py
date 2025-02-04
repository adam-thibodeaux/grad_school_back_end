import gemmi
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem, rdForceFieldHelpers, rdchem, rdmolops, rdMolAlign
# from rdkit.ForceField import rdForceField
import extract_sols as extract
import matplotlib.pyplot as plt
import numpy as np
# import time
# import math
import remove_model
# import extract_sols
import os
import remove_model
import show_dihedral_change
import subprocess
import random
import math
import reorder_pdb
from pathlib import Path
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

# 0,8,9,12,6,11,4,5,

allowable_dihedral_angles_graz = [
    #0 - sulfone side chain flip amide
    (["C30","N3","C29","O5"],[],False),
    #1 - sulfone side chain - carbon part - this is actually not right, idt
    (["C29", "N3","C30","C35"], [],False),
    #2 - sulfone side chain - sulfur part (more medial)
    (["S1","N4","C35","O6"], [],False),
    #3 - methoxy group distal
    (["C24","O9","C22","C23"],[],False),
    #4 - amide link in side chain - shouldn't be necessary (captured in 0)
    (["C1","C29","N3","C30"],[["H5"]],True),
    #5 - sulfone side chain no flip amide
    (["N2","C1","C29","N3"],[],False),
    #6 - sulfone side chain - sulfur part (more distal)
    (["C35","N4","S1","C36"], [],False),
    #7 - macrocycle core distal pyrimidine
    (["C10","C11","C12","C13"], [["H17","H18"]],True),
    #8 - macrocycle core medial pyrimidine
    (["C16","C15","C14","C13"], [["H21","H22"]],True),
    #9 - macrocycle core super medial pyrimidine
    (["N6","C16","C15","C14"], [["H23","H24"]],True)
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
def main(dihedral_angle_atom_names, split_idx_atom_names, count, perturb_angle, input_struct_path, should_proximity_bond=False, should_sanitize=True,should_use_electron_difference=False):
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
        should_valid_struct = True
        # angle = orig_dihedral_angles[i]
        angle_indices = dihedral_angle_atom_indices[i]
        # print(angle)
        # atom_indices = [atom_idx_dict[x] for x in allowable_dihedral_angles[0]]
        comp_id_pos = base_mol.AddConformer(base_mol.GetConformer(0), assignId=True)
        # comp_id_neg = base_mol.AddConformer(base_mol.GetConformer(0), assignId=True)
        edit_mol = rdchem.RWMol(base_mol)
       
        phaser_mtz = input_struct_path.split("/")
        name = phaser_mtz.pop()
        # if "refine" in input_struct_path:
            # phaser_mtz.append(name.replace(".pdb",".mtz"))
        # else:
        phaser_mtz.append(name.split("_out")[0]+"_out.1.mtz")
        phaser_mtz = "/".join(phaser_mtz)
        try:
            # setDihedralForOneBranch(edit_mol, edit_mol.GetConformer(comp_id_pos), atom_indices, perturb_angle)
            rdMolTransforms.SetDihedralDeg(base_mol.GetConformer(comp_id_pos),*angle_indices, perturb_angle)
            # conf = base_mol.GetConformer(comp_id_pos)
            # for bond in base_mol.GetBonds():
            #         start = bond.GetBeginAtomIdx()
            #         end = bond.GetEndAtomIdx()
            #         start_coords = np.array(conf.GetAtomPosition(start))
            #         end_coords = np.array(conf.GetAtomPosition(end))
            #         dist = np.linalg.norm(end_coords-start_coords)
            #         print(dist)
            #         if dist > 3 and should_valid_struct:
            #             print(f'suspicious bond distance {dist} angstroms for: paritaprevir_torsion_angle_perturb_{count[0]}.pdb')
            #             print('skipping this one')
            #             return
            rdForceFieldHelpers.MMFFOptimizeMolecule(base_mol)
            Chem.MolToPDBFile(base_mol, f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count[0]}.pdb", confId=comp_id_pos)
            if should_use_electron_difference:

                return calc_new_fit(f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count[0]}.pdb",phaser_mtz,count, input_struct_path)
            else:  
                # print(count)
                # print([count[0],[0 for i in range(count[0])]])
                return [count[0],0]
            
        except Exception as e:
            print(e)
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
                    if dist > 3 and should_valid_struct:
                        print(f'suspicious bond distance {dist} angstroms for: paritaprevir_torsion_angle_perturb_{count[0]}.pdb')
                        print('skipping this one')
                        return
                # Chem.SanitizeMol(edit_mol, )
            # edit_mol = AllChem.AssignBondOrdersFromTemplate(reference_mol, edit_mol)
            

            # edit_mol.GetConformer(comp_id_pos).SetAtomPosition(bi,ba_coords)
            # edit_mol.GetConformer(comp_id_pos).SetAtomPosition(ei,ea_coords)
            # print(ba_coords.x)
            # print(ea_coords.x)
            Chem.SanitizeMol(edit_mol)
            rdForceFieldHelpers.MMFFOptimizeMolecule(edit_mol)
            
            Chem.MolToPDBFile(edit_mol, f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count[0]}.pdb", confId=comp_id_pos)

            return calc_new_fit(f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count[0]}.pdb", phaser_mtz, count, input_struct_path)
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

def calc_new_fit(pdb_path, mtz_path, count, orig_pdb):
    this_mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
    orig_mol = Chem.MolFromPDBFile(orig_pdb, removeHs=True)

    mtz = gemmi.read_mtz_file(mtz_path)
    size = mtz.get_size_for_hkl()
    # grid = mtz.get_f_phi_on_grid('DELFWT', "PHDELWT", size)
    # # complex_map = np.fft.ifftn(grid.array.conj())
    # # scale_factor = complex_map.size / grid.unit_cell.volume
    # # real_map = np.real(complex_map) * scale_factor
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = mtz.transform_f_phi_to_map("DELFWT", "PHDELWT", size)
    i = 0
    intensities = []
    for _ in this_mol.GetAtoms():
        position = this_mol.GetConformer(0).GetAtomPosition(i)
        orig_position = orig_mol.GetConformer(0).GetAtomPosition(i)
        if abs(orig_position.x - position.x) < 0.01 and abs(orig_position.y - position.y) < 0.01 and abs(orig_position.z - position.z) < 0.01:
            i+= 1
            continue
        intensity = ccp4.grid.interpolate_value(gemmi.Position(position.x,position.y,position.z))
        if intensity < 0:
            intensities.append(abs(intensity))
        i+= 1
    # position = this_mol.GetConformer(cid).GetAtomPosition(atom_idx_input)
    #                 # print(position.x - base_position.x, position.y-base_position.y, position.z-base_position.z)
    # intensity = ccp4.grid.interpolate_value(gemmi.Position(position.x,position.y,position.z))
    # print(sum(intensities))
    # print(intensities)
    return [count[0], sum(intensities)]
    


def shell_source(script):
    """Sometime you want to emulate the action of "source" in bash,
    settings some environment variables. Here is a way to do it."""
    pipe = subprocess.Popen(f"/bin/bash {script}", stdout=subprocess.PIPE, shell=True)
    output = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in output.splitlines()))
    os.environ.update(env)

def generate_conformers(count, step, dihedral_angle_atom_names, split_idx_atom_names, input_struct_path, should_proximity_bond=False):

    intensities = []
    print(dihedral_angle_atom_names)
    # return
    iter_start = 0
    iter_end = 360
    for amide in amide_bonds:
        if dihedral_angle_atom_names[1] in amide and dihedral_angle_atom_names[2] in amide:
            print('bond is an amide, enforcing trans planar geometry')
            print(dihedral_angle_atom_names)
            iter_start = 170
            iter_end = 191
            step = 10
            mol = Chem.MolFromPDBFile(input_struct_path, removeHs=False)
            atom_idx_dict = {}
            
            for atom in mol.GetAtoms():
                atom_idx_dict[atom.GetMonomerInfo().GetName().strip()] = atom.GetIdx()
                atom_idx_dict[str(atom.GetIdx())] = atom.GetMonomerInfo().GetName().strip()

            for neighbor in mol.GetAtomWithIdx(atom_idx_dict[dihedral_angle_atom_names[1]]).GetNeighbors():
                if "C" in dihedral_angle_atom_names[1]:
                    if neighbor.GetAtomicNum() == 8:
                        dihedral_angle_atom_names[0] = atom_idx_dict[str(neighbor.GetIdx())]
                elif "N" in dihedral_angle_atom_names[1]:
                    if neighbor.GetAtomicNum() == 1:
                        dihedral_angle_atom_names[0] = atom_idx_dict[str(neighbor.GetIdx())]
                        
            for neighbor in mol.GetAtomWithIdx(atom_idx_dict[dihedral_angle_atom_names[2]]).GetNeighbors():
                if "C" in dihedral_angle_atom_names[2]:
                    if neighbor.GetAtomicNum() == 8:
                        dihedral_angle_atom_names[3] = atom_idx_dict[str(neighbor.GetIdx())]
                elif "N" in dihedral_angle_atom_names[2]:
                    if neighbor.GetAtomicNum() == 1:
                        dihedral_angle_atom_names[3] = atom_idx_dict[str(neighbor.GetIdx())]
            print(dihedral_angle_atom_names)
            break
            

    for run_num in range(iter_start,iter_end,step):
        count[0] += 1
        try:
            # print(count[0])
            new_intense = main(dihedral_angle_atom_names=dihedral_angle_atom_names, split_idx_atom_names=split_idx_atom_names, count=count, perturb_angle=run_num, input_struct_path=input_struct_path, should_proximity_bond=should_proximity_bond)
            if new_intense:
                intensities.append(new_intense)

            # Chem.MolToPDBFile(mol, f"{output_conf_dir}paritaprevir_torsion_angle_perturb_{count}.pdb", confId=1)
            continue
        except Exception as e:
            print(e)
            continue
    intensities = list(filter(lambda x : x[1] != 0, intensities))
    sorted_intensities = sorted(intensities, key= lambda x : abs(x[1] - 0))
    count[1] += sorted_intensities
    count[1] = [x for x in count[1] if len(x) > 1]
    count[1] = sorted(count[1], key = lambda x: abs(x[1]-0))
    
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
    data.sort_values(["llg","tfz"], ascending=[False,False], inplace=True)
    print(data.to_string())
    set_input_path_list = set()
    list_input_path_list = []
    i = 0
    while len(set_input_path_list) < num_to_include and i < len(data):
        prefix = data.iloc[i]["name"].split("_out")[0]
        if prefix not in set_input_path_list:
            set_input_path_list.add(prefix)
            list_input_path_list.append(data.iloc[i]["name"].split(".")[0])
        i += 1
   
    # return [f{file_prefix}]
    def transform_file_name(name):
        file_name = input_dir
        file_name += name.split("_out")[0]
        file_name += "/"
        file_name += name
        file_name += ".1.pdb"
        return file_name
    return [transform_file_name(x) for x in list_input_path_list]
    ret_list = [input_dir + data.iloc[0]["name"].split("_out")[0] +"/"+ data.iloc[0]["name"]+".pdb",
            input_dir + data.iloc[1]["name"].split("_out")[0] +"/"+ data.iloc[1]["name"]+".pdb",
            input_dir + data.iloc[2]["name"].split("_out")[0] +"/"+ data.iloc[2]["name"]+".pdb"]
    return ret_list[:num_to_include]
def transform_name(input_list, input_mol):
    output_list = []
    for input_item in input_list:
        output_list.append(input_mol.GetAtomWithIdx(input_item).GetMonomerInfo().GetName().strip())
    return output_list
def check_for_overlap(mol, threshold = 2.0):
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        disallowed_idxs = [atom.GetIdx()]
        atom_coords = conf.GetAtomPosition(atom.GetIdx())
        for bond in atom.GetBonds():
            disallowed_idxs.append(bond.GetOtherAtomIdx(atom.GetIdx()))
        for inner_atom in mol.GetAtoms():
            if inner_atom.GetIdx() in disallowed_idxs:
                continue
            inner_atom_coords = conf.GetAtomPosition(inner_atom.GetIdx())
            distance = Chem.rdGeometry.Point3D.Distance(atom_coords, inner_atom_coords)
            if distance < threshold:
                return True
    return False


def rank_dihedrals(input_pdb_path_this):
    mol = Chem.MolFromPDBFile(input_pdb_path_this, removeHs=True)
    dihed_smarts_query = '[!$(*#*)&!D1]~[!$(*#*)&!D1]'
    dihed_smarts_mol = Chem.MolFromSmarts(dihed_smarts_query)
    matches = mol.GetSubstructMatches(dihed_smarts_mol)
    dihed_list = []
    for match in matches:
        idx2 = match[0]
        idx3 = match[1]
        bond = mol.GetBondBetweenAtoms(idx2, idx3)
        jAtom = mol.GetAtomWithIdx(idx2)
        kAtom = mol.GetAtomWithIdx(idx3)
        for b1 in jAtom.GetBonds():
            if (b1.GetIdx() == bond.GetIdx()):
                continue
            idx1 = b1.GetOtherAtomIdx(idx2)
            for b2 in kAtom.GetBonds():
                if ((b2.GetIdx() == bond.GetIdx())
                or (b2.GetIdx() == b1.GetIdx())):
                    continue
                idx4 = b2.GetOtherAtomIdx(idx3)
            # skip 3-membered rings
                if (idx4 == idx1):
                    continue
                dihed_list.append((idx1, idx2, idx3, idx4))
    
    allowed_dihed_angles = []
    for angle in dihed_list:
        try:
            orig_angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(),angle[0],angle[1],angle[2],angle[3])
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(),angle[0],angle[1],angle[2],angle[3], orig_angle)
        except:
            continue
        allowed_dihed_angles.append(angle)
        # print([transform_name(x, mol) for x in allowed_dihed_angles])
    mol = Chem.MolFromPDBFile(input_pdb_path_this)
    for i in range(len(allowed_dihed_angles)):
        edit_mol = Chem.MolFromPDBFile(input_pdb_path_this)
        dihed = allowed_dihed_angles[i]
        rmsd_list = []
        for j in range(0,360,10):
            rdMolTransforms.SetDihedralDeg(edit_mol.GetConformer(),dihed[0],dihed[1],dihed[2],dihed[3],j)
            if check_for_overlap(edit_mol):    
                # print('faulty conformer found with invalid bonding, adding penalty to RMSD calculation')
                rmsd_list.append(-10.0)
                continue
            # print(f"after perturbation, edit_mol has {len(edit_mol.GetBonds())} bonds")

            rmsd = rdMolAlign.GetBestRMS(mol,edit_mol)
            if rmsd > 100000:
                print('invalid torsion angle found, applying penalty and skipping')
                rmsd_list.append(-10.0)
                continue
            rmsd_list.append(rmsd)
        allowed_dihed_angles[i] = (dihed, transform_name(dihed,mol), rmsd_list, sum(rmsd_list)/len(rmsd_list), False)
    # print(f"x[-1]\n" for x in allowed_dihed_angles])
    allowed_dihed_angles = sorted(allowed_dihed_angles, key=lambda x : x[-2],reverse=True)
    non_dupes = []
    for final in allowed_dihed_angles:
        should_append = True
        for non_dupe in non_dupes:
            if abs(final[-2] - non_dupe[-2]) < 0.005:
                print('duplicate dihedral found')
                print(non_dupe[1],non_dupe[-2])
                print(final[1],final[-2])
                print(abs(final[-2] - non_dupe[-2]))
                should_append = False
                break
        if should_append:
            non_dupes.append(final)
    for non_dupe in non_dupes:
        print(non_dupe[1])
        print(non_dupe[-2])
    non_dupes_amide_ordered = []
    skip_idx = []
    for j in range(len(non_dupes)):
        non_dupe = non_dupes[j]
        print(non_dupe[1])
        for amide in amide_bonds:
            if non_dupe[1][1] in amide and non_dupe[1][2] in amide:
                print('bond is an amide, ranking first')
                # print(non_dupe)
                
                non_dupes_amide_ordered.append(non_dupe)
                skip_idx.append(j)
    # print(amide_bonds)
    # print(non_dupes_amide_ordered)
    for j in range(len(non_dupes)):
        if j in skip_idx:
            continue
        non_dupes_amide_ordered.append(non_dupes[j])

    return non_dupes_amide_ordered

def get_custom_diheds(input_struct_path_this_one):

    mol = Chem.MolFromPDBFile(input_struct_path_this_one, removeHs=False)
    atom_idx_dict = {}
    for atom in mol.GetAtoms():
        atom_idx_dict[atom.GetMonomerInfo().GetName()] = atom.GetIdx()
    ring_info = mol.GetRingInfo()
    ring_systems = []
    for ring in ring_info.AtomRings():
        ring_atoms = set(ring)
        this_ring_systems = []
        for ring_system in ring_systems:
            num_overlap = len(ring_atoms.intersection(ring_system))
            if num_overlap > 1:
                ring_atoms = ring_atoms.union(ring_system)
            else:
                this_ring_systems.append(ring_system)
        this_ring_systems.append(ring_atoms)
        ring_systems = this_ring_systems
    print(ring_systems)
    largest_system = list(sorted(ring_systems, key= lambda x: len(x), reverse=True))[0]
    ret_list = []
    for atom_idx in largest_system:
        is_three_member_ring = False
        target_atom = mol.GetAtomWithIdx(atom_idx)
        if target_atom.GetIsAromatic() or (target_atom.GetHybridization() == Chem.HybridizationType.SP2 and target_atom.GetTotalNumHs(includeNeighbors=True) > 0):
            print('aromatic or sp2 macrocycle atom found, skipping custom perturbation')
            continue
        #fetch non-terminal neighbors: if > 2 cant do it
        neighbors = target_atom.GetNeighbors()
        sub_neighbors = []
        terminal_neighbors = []
        sub_sub_neighbor = None
        for neighbor in neighbors:
            this_sub_neighbors = neighbor.GetNeighbors()
            if len(this_sub_neighbors) > 1:
                sub_neighbors.append(neighbor.GetMonomerInfo().GetName().strip())
            else:
                terminal_neighbors.append(neighbor.GetMonomerInfo().GetName().strip())
            if sub_sub_neighbor == None:
                for sub_neighbor in this_sub_neighbors:
                    if sub_neighbor.GetMonomerInfo().GetName().strip() == target_atom.GetMonomerInfo().GetName().strip():
                        continue
                    this_sub_sub_neighbors = sub_neighbor.GetNeighbors()
                    for this_sub_sub_neighbor in this_sub_sub_neighbors:
                        if sub_sub_neighbor != None: 
                            continue
                        if this_sub_sub_neighbor.GetMonomerInfo().GetName().strip() == target_atom.GetMonomerInfo().GetName().strip():
                            continue
                        sub_sub_sub_neighbors = this_sub_sub_neighbor.GetNeighbors()
                        if len(sub_sub_sub_neighbors) > 1:
                            sub_sub_neighbor = (this_sub_sub_neighbor.GetMonomerInfo().GetName().strip(), sub_neighbor.GetMonomerInfo().GetName().strip())
        if len(terminal_neighbors) == 0 or len(sub_neighbors) > 2 or sub_sub_neighbor == None:
            continue
        if (len(sub_neighbors) + len(terminal_neighbors))  == 3 and len(terminal_neighbors) > 0 and terminal_neighbors[0][0] == "H" and target_atom.GetMonomerInfo().GetName().strip()[0] == "C":
            continue
        # print(f"terminal neighbors for {target_atom.GetMonomerInfo().GetName()}")
        # print(terminal_neighbors)
        # print(f"non-terminal neighbors for {target_atom.GetMonomerInfo().GetName()}")
        # print(sub_neighbors)
        # print(f"final neighbor for {target_atom.GetMonomerInfo().GetName().strip()}")
        # print(sub_sub_neighbor)
        out_neighbor = None
        for sub_neighbor in sub_neighbors:
            if sub_neighbor != sub_sub_neighbor[1] and sub_neighbor != sub_sub_neighbor[0]:
                out_neighbor = sub_neighbor
                break
        if not is_three_member_ring:
            ret_list.append(([],[sub_sub_neighbor[0],sub_sub_neighbor[1],target_atom.GetMonomerInfo().GetName().strip(), out_neighbor],terminal_neighbors, 0,  True))
    print(ret_list)
    return ret_list
def find_amides(mol):
    print('finding amides')
    amide_smiles = "C(=O)-N"
    template_smiles = "Cc1cnc(cn1)C(=O)N[C@H]2CCCCC/C=C\\[C@@H]3C[C@]3(NC(=O)[C@@H]4C[C@H](CN4C2=O)Oc5c6ccccc6c7ccccc7n5)C(=O)NS(=O)(=O)C8CC8"
    amide_query = Chem.MolFromSmarts(amide_smiles)
    
    template = Chem.MolFromSmiles(template_smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    matches = mol.GetSubstructMatches(amide_query)
    atom_idx_dict = {}
    for atom in mol.GetAtoms():
        atom_idx_dict[atom.GetIdx()] = atom.GetMonomerInfo().GetName().strip()
    matches = [[atom_idx_dict[match[0]],atom_idx_dict[match[1]],atom_idx_dict[match[2]]] for match in matches]
    for match in matches:
        print(match)
    return matches
    
            




input_struct_path_list = []
# perturb_path = [x for x in range(len(allowable_dihedral_angles_sim)) if allowable_dihedral_angles_sim[x] != "REFINE"]
# temp = [*perturb_path]
# perturb_path = [4]
# perturb_path += temp


# [0,8,3,7] give a somewhat nice structure 
perturb_path = [0,8,7,9,3]
perturb_path =[1]

# for i in range(100):
#     num = random.randrange(len(sub_allowed))
#     perturb_path.append(sub_allowed[num])
    
input_struct_path = "/Users/adam/Downloads/inputs_for_molec_replac/PAR_BETA_CONFORGE_TRIAL_1/paritaprevir_beta_conforge_78.pdb"
input_mtz = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta_20A.mtz"
monomer_library = "/Users/adam/Downloads/outputs_from_molec_replac/phenix_refine/eLBOW_65/elbow.UNK.001.cif"
output_path = "/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CUSTOM_CONF_TRIAL_5_20A"
Path(output_path).mkdir(exist_ok=True, parents=True)
with open(input_struct_path,"r") as f:
    string = f.read()
    input_struct_name = input_struct_path.split("/")[-1]
    with open(f"{output_path}/{input_struct_name}","w") as new_f:
        new_f.write(string)
        new_f.close()
    f.close()


# perturb_path = get_custom_diheds(input_struct_path_this_one=input_struct_path)
# perturb_path.append("REFINE")
amide_bonds = find_amides(Chem.MolFromPDBFile(input_struct_path))
perturb_path = rank_dihedrals(input_pdb_path_this=input_struct_path)
perturb_path = list(filter(lambda x: x[-2] > -3, perturb_path))


# perturb_path = [*perturb_path[1:]]
perturb_path.append("REFINE")
print(perturb_path)
# print(len(perturb_path))
for i in range(len(perturb_path)):
    
    # AllChem.AssignBondOrdersFromTemplate(tempalte, mol)

    print(input_struct_path_list)
    # print(perturb_path[i][-2])
    print(perturb_path[i])
    #skipping dihedral angles that cause lots of overlap
    if perturb_path[i] != "REFINE" and perturb_path[i][-2] < -3: 
        continue
    output_conf_dir = f"{output_path}/ROUND_{i}/"
    count = [0,[]]
    should_skip_conf_gen = False
    should_use_electron_difference = False
    num_to_include = 5
    if os.path.isdir(output_conf_dir) and len(os.listdir(output_conf_dir)) > 2:
        print(f"round {i} already complete, skipping conformer_generation")
        should_skip_conf_gen = True
    else:

        os.makedirs(output_conf_dir, exist_ok=True)
        os.makedirs(output_conf_dir + "PHASER", exist_ok=True)
    if perturb_path[i-1] == "REFINE":
        should_use_electron_difference = False

    if perturb_path[i] == "REFINE":
        
        for k in range(len(input_struct_path_list)):
            input_struct_path_ind = input_struct_path_list[k]
            
            
            fixed_input_struct_path = output_conf_dir.split("ROUND")[0]
            fixed_input_struct_path += f"ROUND_{i-1}/PHASER/"
            input_file_name = input_struct_path_ind.split("/")[-1].split("_out")[0]
            fixed_input_struct_path += input_file_name
            fixed_input_struct_path += "/"
            fixed_input_struct_path += input_file_name
            fixed_input_struct_path += "_out.1.pdb"

            source_str = "/usr/local/phenix-1.21.2-5419/phenix_env.sh"
            prefix = input_struct_path_ind.split("/")[-1].split(".pdb")[0]
            # mol_ind = Chem.MolFromPDBFile(fixed_input_struct_path, removeHs=True, sanitize=False, proximityBonding=False)
            # mol_ind = Chem.RemoveAllHs(mol_ind)
            # Chem.MolToPDBFile(mol_ind, f"{output_conf_dir}/{prefix}.pdb")
            # fixed_input_struct_path = f"{output_conf_dir}/{prefix}.pdb"
            # print(f"{prefix}_refine_001.mtz")
            if f"{prefix}_refine_001.pdb" in os.listdir(output_conf_dir):
                print('refinement already done for this structure, skipping')
                continue
            print(f'starting initial refinement for {prefix}')
            param_str = ""
            with open("/Users/adam/Downloads/inputs_for_molec_replac/new_phenix_params_v2","r") as f:
                param_str = f.read()
                f.close()
            param_str = param_str.replace("REPLACE_PDB", fixed_input_struct_path)
            param_str = param_str.replace("REPLACE_MTZ", input_mtz)
            param_str = param_str.replace("REPLACE_RESTRAINT", monomer_library)
            param_str = param_str.replace("REPLACE_SPACEGROUP", "P 21 21 21")
            param_str = param_str.replace("REPLACE_UNIT_CELL", "10.56 12.32 31.73 90.00 90.00 90.00")
            param_str = param_str.replace("REPLACE_RFLAG_FRACTION", "0.05")
            with open(f"{output_conf_dir}/phenix_params_{k}","w") as f:
                f.write(param_str)
                f.close()
            #run refine on each of these input_structs
            shell_source(source_str)
            proc = subprocess.run([f"cd {output_conf_dir} && source {source_str} && phenix.refine ./phenix_params_{k} overwrite=true refinement.input.symmetry_safety_check=warning" ], capture_output=True, shell=True)
            print(proc.stderr.decode().strip())
            print(proc.stdout.decode().strip())
            print(f'finished initial refinement for {prefix}')
            # print(f"adding hydrogens to refined structure")
            # to_remove_atoms = ["H13","H14"]
            # fixed_dihedral = ["H13","C9","C10","H14"]
            # mol_ref = Chem.MolFromPDBFile(f"{output_conf_dir}/{prefix}_refine_001.pdb", removeHs=False)
            # coord_dict = {}
            # conf = mol_ref.GetConformer()
            # atom_idx_dict = {}
            # for atom in mol_ref.GetAtoms():
            #     atom_name = atom.GetMonomerInfo().GetName().strip()
            #     atom_idx = atom.GetIdx()
            #     if atom_name not in to_remove_atoms:
            #         coord_dict[atom.GetIdx()] = conf.GetAtomPosition(atom_idx)
            #     atom_idx_dict[atom_name] = atom_idx
            # mp =  rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_ref)
            # ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_ref, mp)
            # for coord in coord_dict.items():
            #     ff.MMFFAddPositionConstraint(coord[0], 0, 1.e4)
            # ff.MMFFAddTorsionConstraint(atom_idx_dict[fixed_dihedral[0]],atom_idx_dict[fixed_dihedral[1]],atom_idx_dict[fixed_dihedral[2]],atom_idx_dict[fixed_dihedral[3]],False, 0.0,0.0,1.e4)
            # ff.Minimize(maxIts=10000)
            # Chem.MolToPDBFile(mol_ref,f"{output_conf_dir}/{prefix}_refined_alkene_fix.pdb")
            # reorder_pdb.reorder_atoms_v2(f"{output_conf_dir}/{prefix}_refined_with_H.pdb", input_struct_path, f"{output_conf_dir}/{prefix}_refined_with_H.pdb")
            # with open("/Users/adam/Downloads/inputs_for_molec_replac/new_phenix_params_electron_v2","r") as f:
            #     param_str = f.read()
            #     f.close()
            # param_str = param_str.replace("REPLACE_PDB", f"{output_conf_dir}/{prefix}_refined_alkene_fix.pdb")
            # param_str = param_str.replace("REPLACE_MTZ", input_mtz)
            # param_str = param_str.replace("REPLACE_RESTRAINT", monomer_library)
            # param_str = param_str.replace("REPLACE_SPACEGROUP", "P 21 21 21")
            # param_str = param_str.replace("REPLACE_UNIT_CELL", "5.090 15.61 50.78 90.00 90.00 90.00")
            # param_str = param_str.replace("REPLACE_RFLAG_FRACTION", "0.05")
            # with open(f"{output_conf_dir}/phenix_params_electron_{k}","w") as f:
            #     f.write(param_str)
            #     f.close()
            # shell_source(source_str)
            # print("starting second round of refinement")
            # proc = subprocess.run([f"cd {output_conf_dir} && source {source_str} && phenix.refine ./phenix_params_electron_{k} overwrite=true refinement.input.symmetry_safety_check=warning" ], capture_output=True, shell=True)
            # print(proc.stderr.decode().strip())
            # print(proc.stdout.decode().strip())
            # print(f'finished second round of refinement for {prefix}')
        input_struct_path_list = []
        for f in os.listdir(output_conf_dir):
            if ".pdb" in f:
                input_struct_path_list.append(f"{output_conf_dir}{f}")
        continue
    dihedral_angle = perturb_path[i][1]
    
    should_use_custom_dihed = perturb_path[i][-1]
    split_atom_names =[]
    if should_use_custom_dihed:
        split_atom_names = [perturb_path[i][-3]]

    # split_atom_names = allowable_dihedral_angles_graz[perturb_path[i]][1]
    print(dihedral_angle)
    # if i == 19:
        # input_struct_path_list = ["/Users/adam/Downloads/outputs_from_molec_replac/PAC_CUSTOM_CONF_TRIAL_1_20A/ROUND_18/paritaprevir_torsion_angle_perturb_6_out.1_refine_001.pdb"]
        # should_use_electron_difference = False
    if i > 0:
        for input_struct_path_individ in input_struct_path_list:
            if should_skip_conf_gen: continue
            input_struct_path = input_struct_path_individ
            print(f"starting conformer generation for {input_struct_path}")
            # if i > 1 and allowable_dihedral_angles_pac[perturb_path[i-1]] == "REFINE":
                # count = generate_conformers(count=count, step=20, dihedral_angle_atom_names=dihedral_angle,split_idx_atom_names=split_atom_names, input_struct_path=input_struct_path, should_proximity_bond=True)
            if i <= 1 and not should_use_custom_dihed:
                count = generate_conformers(count=count, step=20, dihedral_angle_atom_names=dihedral_angle,split_idx_atom_names=split_atom_names, input_struct_path=input_struct_path, should_proximity_bond=True)
            elif should_use_custom_dihed:
                count = show_dihedral_change.main(input_struct = input_struct_path_individ, atom_names = dihedral_angle,adjacent_atom_names=split_atom_names[0],output_path=output_conf_dir, count = count, num_points=20, mtz_path=input_mtz)
            else:
                count = generate_conformers(count=count, step=20, dihedral_angle_atom_names=dihedral_angle,split_idx_atom_names=split_atom_names, input_struct_path=input_struct_path, should_proximity_bond=True)
        timers = []
        sorted_output_conf = []
        for j in os.listdir(output_conf_dir):
            
            if should_use_electron_difference:
                allowed_conformers = [str(x[0]) for x in count[1][:num_to_include]]
                # print(allowed_conformers)
                if (j.split("_")[-1].split(".pdb")[0]) in allowed_conformers:
                    sorted_output_conf.append(j)
            else:
                sorted_output_conf.append(j)
        sorted_output_conf = sorted(sorted_output_conf)
        if should_use_electron_difference:
            print("expected top 5 according to electron difference map")
            print(count[1][:5])
            if len(count[1]) > 0:
                with open(f"{output_conf_dir}electron_difference_summary.txt","w") as file:
                    file.write(str(count))
                    file.close()
        for k in range(len(sorted_output_conf)):
            j = sorted_output_conf[k]
            phaser_input = output_conf_dir + j
            if ".pdb" not in phaser_input: continue
            phaser_name = phaser_input.split("/")[-1].replace(".pdb","")
            if os.path.isdir(output_conf_dir + "PHASER/" + phaser_name):
                print(f'PHASER already done for {phaser_name}, skipping, round {i}')
                continue

            print('starting PHASER for ' + phaser_name)
            timer = remove_model.main_tfd(input_pdb=phaser_input,output_dir=output_conf_dir+"PHASER/"+phaser_name+"/", out_dir= phaser_name, input_mtz_tfd=input_mtz)
            timers.append(timer)
            print(f'finished PHASER for {phaser_name} in {timer} seconds')
            print(f'approximate time remaining for this round {i} is: {(sum(timers)/len(timers)*(len(sorted_output_conf) - k - 1))/60} minutes')
        input_struct_path_list = extract_llg_and_tfz_and_update_init_path(output_conf_dir + "PHASER/", num_to_include=num_to_include)

    else:
        print(f"starting initial conformer generation")
        if should_use_custom_dihed:
            count = show_dihedral_change.main(input_struct = input_struct_path, atom_names = dihedral_angle,adjacent_atom_names=split_atom_names[0],output_path=output_conf_dir, count = count, num_points=20, mtz_path=input_mtz)
        else:
            generate_conformers(count=count, step=20, dihedral_angle_atom_names=dihedral_angle, split_idx_atom_names=split_atom_names, input_struct_path=input_struct_path, should_proximity_bond=True)
        input_struct_path_list = output_conf_dir
        timers = []
        sorted_output_conf = []
        for j in os.listdir(output_conf_dir):
            sorted_output_conf.append(j)
        sorted_output_conf = sorted(sorted_output_conf)
        for k in range(len(sorted_output_conf)):
            j = sorted_output_conf[k]
            phaser_input = output_conf_dir + j
            if ".pdb" not in phaser_input: continue
            phaser_name = phaser_input.split("/")[-1].replace(".pdb","")
            if os.path.isdir(output_conf_dir + "PHASER/" + phaser_name):
                print(f'PHASER already done for {phaser_name}, skipping, round {i}')
                continue
            print('starting PHASER for ' + phaser_name)
            timer = remove_model.main_tfd(input_pdb=phaser_input,output_dir=output_conf_dir+"PHASER/"+phaser_name+"/", out_dir= phaser_name, input_mtz_tfd=input_mtz)
            timers.append(timer)
            print(f'finished PHASER for {phaser_name} in {timer} seconds')
            print(f'approximate time remaining for this round {i} is: {((sum(timers)/len(timers)) * (len(sorted_output_conf) - k - 1))/60} minutes')
        input_struct_path_list = extract_llg_and_tfz_and_update_init_path(output_conf_dir + "PHASER/", num_to_include=num_to_include)
        if should_use_electron_difference:
            print("expected top 5 according to electron difference map")
            print(count[1][:5])
            if len(count[1]) > 0:
                with open(f"{output_conf_dir}electron_difference_summary.txt","w") as file:
                    file.write(str(count))
                    file.close()

        


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
# #     if ".pdb" not in f.split(".")[-1]: continue
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

        

