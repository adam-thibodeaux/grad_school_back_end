from Bio.PDB import PDBParser, Superimposer, Selection
import pymatgen.command_line
from rdkit.Chem import rdMolAlign, TorsionFingerprints, AllChem, rdMolTransforms, rdForceFieldHelpers
from rdkit import Chem
from rdkit.Geometry import Point3D
import functools

import numpy as np
import math
import scipy.spatial.transform as transform 

# def rotation_matrix(axis, theta):
#     """
#     Return the rotation matrix associated with counterclockwise rotation about
#     the given axis by theta radians.
#     """
#     axis = np.asarray(axis)
#     axis = axis / math.sqrt(np.dot(axis, axis))
#     a = math.cos(theta / 2.0)
#     b, c, d = -axis * math.sin(theta / 2.0)
#     aa, bb, cc, dd = a * a, b * b, c * c, d * d
#     bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
#     return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
#                      [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
#                      [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
def calc_rmsd_pdb(input_struct_dir, compare_struct_dir):
    parser = PDBParser()
    superimposer = Superimposer()
    input_struct = parser.get_structure(id="INP", file=input_struct_dir)
    compare_struct = parser.get_structure(id="COM", file=compare_struct_dir)
    atoms_1 = Selection.unfold_entities(input_struct,"A")
    atoms_2 = Selection.unfold_entities(compare_struct,"A")
    superimposer.set_atoms(atoms_1, atoms_2)
    superimposer.apply(compare_struct.get_atoms())
    return superimposer.rms
    
def calc_rmsd(input_struct_path, compare_struct_path):
    input_mol = Chem.MolFromPDBFile(input_struct_path)
    compare_mol = Chem.MolFromPDBFile(compare_struct_path)
    # input_mol = AllChem.AssignBondOrdersFromTemplate(compare_mol, input_mol)
    return rdMolAlign.GetBestRMS(input_mol, compare_mol)


def calc_tfd(input_struct_path, compare_struct_path):
    smiles_str = "C[C@@H]1CNC(C(O)N[C@H]2CCCCCCC[C@@H]3C[C@@]3(C(O)NS(O)(O)C3CC3)NC(O)[C@@H]3C[C@@H](O[C@@H]4NC5CCCC[C@H]5[C@@H]5CCCCC54)CN3[C@H]2O)CN1"
    reorder_mol = Chem.MolFromSmiles(smiles_str)
    # print(len(reorder_mol.GetAtoms()))
    input_mol = Chem.MolFromPDBFile(input_struct_path)
    input_mol = AllChem.AssignBondOrdersFromTemplate(reorder_mol, input_mol)
    input_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(input_mol))
    # print(len(input_mol.GetAtoms()))
    compare_mol = Chem.MolFromPDBFile(compare_struct_path)
    compare_mol = AllChem.AssignBondOrdersFromTemplate(reorder_mol, compare_mol)
    compare_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(compare_mol))
    # input_match = input_mol.GetSubstructMatch(reorder_mol)
    # print(input_match)
    # input_mol = Chem.RenumberAtoms(input_mol, input_match)
    
    # compare_match = input_mol.GetSubstructMatch(reorder_mol)
    # compare_mol = Chem.RenumberAtoms(compare_mol, compare_match)
    # print(compare_match)

    print(Chem.MolToSmiles(input_mol))
    print(len(input_mol.GetBonds()))
    print(Chem.MolToSmiles(compare_mol))
    print(len(compare_mol.GetBonds()))
    # input_mol = AllChem.AssignBondOrdersFromTemplate(compare_mol, input_mol)
    return TorsionFingerprints.GetTFDBetweenMolecules(input_mol, compare_mol)

if __name__ == "__main__":
    # print(calc_rmsd("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta_reordered.pdb",
    #                 "/Users/adam/Downloads/outputs_from_molec_replac/phenix_refine/Refine_79/parateprevir_refine_079.pdb"
    #                 ))
    # input_mol_no_H = Chem.MolFromSmiles(
    #     "O[C@H]1C[C@@H](CSCCC([O-])=O)O[C@H](O)[C@@H]1O"
    #     )
    # input_mol_with_H = Chem.AddHs(input_mol_no_H)
    # Chem.MolToMolFile(input_mol_with_H, '/Users/adam/Downloads/inputs_for_molec_replac/sugammadex_sub_unit.mol')
    # test_str = """
    # i2run phaser_simple --projectName MR \\
	# --F_SIGF \\
	# 	"project=44190fb09bfb11efa7bd321be340d2e4" \\
	# 	"baseName=F_SIGF-observed_data_19_1.mtz" \\
	# 	"relPath=CCP4_IMPORTED_FILES" \\
	# 	"annotation=Reflections imported from F_SIGF-observed_data_19.mtz by job 358" \\
	# 	"subType=1" \\
	# 	"contentFlag=3" \\
	# --ENSEMBLES \\
	# 	"label=Ensemble_1" \\
	# 	"number=1" \\
	# 	"use=True" \\
	# --SGALT_SELECT ALL \\
	# --RUNSHEETBEND False \\
	# --RUNREFMAC False \\
	# --XYZIN \\
	# 	"project=44190fb09bfb11efa7bd321be340d2e4" \\
	# 	"baseName=XYZIN-coordinates_11_1_1.pdb" \\
	# 	"relPath=CCP4_IMPORTED_FILES" \\
	# 	"annotation=Search imported from XYZIN-coordinates_11_1.pdb by job 358" \\
	# 	"subType=0" \\
	# 	"contentFlag=1" \\
	# --SEARCHSEQUENCEIDENTITY 1.0 \\
	# --XYZIN_FIXED \\
	# 	"project=44190fb09bfb11efa7bd321be340d2e4" \\
	# 	"baseName=PHASER.1.pdb" \\
	# 	"relPath=CCP4_JOBS/job_67" \\
	# 	"annotation=Positioned coordinates for solution 1" \\
	# 	"dbFileId=cc79ce92a2dd11ef91f2321be340d2e4" \\
	# 	"subType=0" \\
	# 	"contentFlag=1" \\
	# --FORM electron \\
	# --PACK_CUTO 100 \\
	# --PACK_QUIC False \\
	# --PACK_COMP False \\
	# --PACK_KEEP_HIGH_TFZ True \\
	# --PEAK_ROTA_SELE all \\
	# --PEAK_TRAN_SELE all \\
	# --PURG_ROTA_ENAB False \\
	# --PURG_TRAN_ENAB False \\
	# --PURG_RNP_ENAB False \\
	# --RESC_ROTA False \\
	# --RESC_TRAN False \\
	# --SEAR_METH full \\
	# --ZSCO_USE False \\
	# --ZSCO_SOLV 4.0 \\
	# --ZSCO_POSS_SOLV 3.0 \\
    # """
    # print(test_str.replace("\n","").replace("\\"," "))
    #a,b,c is point on axis
    #u,v,w is axis vector direction
    #x,y,z is coordinates
    def find_normal(points):
        p0, p1, p2 = points
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        x2, y2, z2 = p2

        ux, uy, uz = u = [x1-x0, y1-y0, z1-z0] #first vector
        vx, vy, vz = v = [x2-x0, y2-y0, z2-z0] #sec vector

        u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx] #cross product

        point  = np.array(p1)
        normal = np.array(u_cross_v)
        return normal
    def rot_point_from_formula(a,  b,  c, u,  v,  w, x,  y,  z, theta):
        
        u2 = u*u
        v2 = v*v
        w2 = w*w
        cosT = np.cos(theta)
        oneMinusCosT = 1 - cosT
        sinT = np.sin(theta)


        p = []
        p.append(((a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * oneMinusCosT
                + x*cosT
                + (-c*v + b*w - w*y + v*z)*sinT))

        p.append(((b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * oneMinusCosT
                + y*cosT
                + (c*u - a*w + w*x - u*z)*sinT))

        p.append(((c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * oneMinusCosT
                + z*cosT
                + (-b*u + a*v - v*x + u*y)*sinT))

        return p

    def normalize(v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm
    center = [0,15,1]
    axis_direction = normalize([0,0,1])

    angle = 0
    
    input_struct = "/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_2/sugammadex_subunit_conforge_4.pdb"
    # core_struct = "/Users/adam/Downloads/inputs_for_molec_replac/sugammadex_core.mol"
    mol = Chem.MolFromPDBFile(input_struct, removeHs=False)
    
    core_mol = Chem.MolFromSmarts("S-C-C-1O-C-C-C-C-1")
    core_mol = Chem.RemoveAllHs(core_mol)
    print(Chem.MolToSmiles(core_mol))
    centroid = rdMolTransforms.ComputeCentroid(mol.GetConformer())
    vol = 20
    center = [centroid.x, centroid.y, centroid.z + vol]
    edit_mol = Chem.RWMol(mol)
    atom_idx_matches_start =[]
    atom_idx_matches_end =[]
    
    prev_mol = edit_mol
    bond_start = "O5"
    bond_end = "C2"
    plane_atoms = ["O4","C9","C2"]
    plane_atom_coords = []
    plane_norm_vec = []
    ring_atoms = ["O4","C9","C2","C1","C3","C8"]
    ring_atom_coords = []
    ring_center = []
    axis_direction_atoms = ["C3","C9"]
    axis_direction_coords = []

    mols = []

    for atom in mol.GetAtoms():
        if atom.GetMonomerInfo().GetName().strip() in plane_atoms:
            position = mol.GetConformer().GetAtomPosition(atom.GetIdx())

            plane_atom_coords.append([position.x, position.y, position.z])
        if atom.GetMonomerInfo().GetName().strip() in ring_atoms:
            position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            ring_atom_coords.append([position.x, position.y, position.z])
        if atom.GetMonomerInfo().GetName().strip() in axis_direction_atoms:
            position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            axis_direction_coords.append([position.x, position.y, position.z])
    for i in range(3):
        ring_center.append(np.average([ring_atom_coords[0][i],ring_atom_coords[1][i],ring_atom_coords[2][i]]))
    plane_norm_vec = find_normal(plane_atom_coords)
    print(ring_center, plane_norm_vec)
    axis_direction = normalize(np.subtract(axis_direction_coords[0],axis_direction_coords[1]))
    center = np.add(np.array(ring_center), (5.5/np.linalg.norm(plane_norm_vec)*np.array(plane_norm_vec)))
    for i in range(8):
        angle += (2*np.pi/8)
        # rotation = transform.Rotation.from_quat([np.cos(angle/2), center[0]* np.sin(angle/2),center[1]* np.sin(angle/2),center[2]* np.sin(angle/2)])
        
        mol = Chem.MolFromPDBFile(input_struct, removeHs=True)

        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            position = conf.GetAtomPosition(atom.GetIdx())
            new_position = rot_point_from_formula(center[0],center[1],center[2],axis_direction[0],axis_direction[1],axis_direction[2],position.x, position.y, position.z, angle)
            # print(position.x,position.y,position.z)
            
            conf.SetAtomPosition(atom.GetIdx(), Point3D(new_position[0], new_position[1], new_position[2]))
        mols.append(mol)

    combined_mols = functools.reduce(Chem.CombineMols, mols)
    for atom in combined_mols.GetAtoms():
        if atom.GetMonomerInfo().GetName().strip() == bond_start:
            atom_idx_matches_start.append(atom.GetIdx())
        elif atom.GetMonomerInfo().GetName().strip() == bond_end:
            atom_idx_matches_end.append(atom.GetIdx())
    print(atom_idx_matches_end)
    print(atom_idx_matches_start)
    combined_mols = Chem.RWMol(combined_mols)
    match_structs = combined_mols.GetSubstructMatches(core_mol)
   

    

    # mp =  rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
            # ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_ref, mp)
            # for coord in coord_dict.items():
            #     ff.MMFFAddPositionConstraint(coord[0], 0, 1.e4)
            # ff.MMFFAddTorsionConstraint(atom_idx_dict[fixed_dihedral[0]],atom_idx_dict[fixed_dihedral[1]],atom_idx_dict[fixed_dihedral[2]],atom_idx_dict[fixed_dihedral[3]],False, 0.0,0.0,1.e4)
            # ff.Minimize(maxIts=10000)
    
    
    # rdForceFieldHelpers.MMFFOptimizeMolecule(combined_mols, maxIters=1000000)
    
    for i in range(len(atom_idx_matches_start) - 1):
        combined_mols.AddBond(atom_idx_matches_start[i], atom_idx_matches_end[i+1], Chem.rdchem.BondType.SINGLE)
    
    combined_mols.AddBond(atom_idx_matches_start[-1], atom_idx_matches_end[0], Chem.rdchem.BondType.SINGLE)
    combined_mols = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromMolFile("/Users/adam/Downloads/inputs_for_molec_replac/sugammadex.mol"),combined_mols)
    Chem.SanitizeMol(combined_mols)
    combined_mols = Chem.AddHs(combined_mols, addCoords=True)
    conf = combined_mols.GetConformer()
    mp =  rdForceFieldHelpers.MMFFGetMoleculeProperties(combined_mols)
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(combined_mols, mp)
    for match_indices in match_structs:
        print(match_indices)
        for match_idx in match_indices:
            ff.MMFFAddPositionConstraint(match_idx, 0.1, 1.e4)
            # continue
    # rdForceFieldHelpers.MMFFOptimizeMolecule(combined_mols, maxIters=1000000)
    Chem.SanitizeMol(combined_mols)
    ff.Minimize(maxIts=10000)
    Chem.MolToMolFile(combined_mols, "./test.mol")
    
