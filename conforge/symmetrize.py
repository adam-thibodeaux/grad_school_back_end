import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdForceFieldHelpers, AllChem, Atom, AtomMonomerInfo, rdqueries
from rdkit.Geometry import Point3D
import functools
import os

def symmetrize(input_struct, output_path, constraint=0.08, center_offset=4.15, axis_offset=0):
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
    def add_sodium(mol):
        anion_oxygen = "O4"
        subunit_struct = Chem.MolFromSmiles("SCCC(O)O")
        subunit_struct = Chem.RemoveAllHs(subunit_struct)
        # print(Chem.MolToSmiles(subunit_struct))
        mol = Chem.RWMol(mol)
        # matches = mol.GetSubstructMatch(subunit_struct)
        # print("length of matches", len(matches))
        should_add_sodium = False
        oxygen_coords_total = []
        count = 0
        ammonium_mol = Chem.MolFromSmiles("[NH4+]")
        print('done')
        bare_nitrogen = Chem.MolFromSmiles("N")
        bare_nitrogen = Chem.RemoveAllHs(bare_nitrogen)
        mol = functools.reduce(Chem.CombineMols, [ammonium_mol, mol])
        mol = Chem.RWMol(mol)

        for atom in mol.GetAtoms():
            # print(atom.GetMonomerInfo().GetName())
            # print(match_idx)
            # atom = mol.GetAtomWithIdx(match_idx)
        # print(atom.GetMonomerInfo().GetName().strip())
        # print(atom.GetExplicitValence())
            if not atom.GetMonomerInfo():
                continue
            if atom.GetMonomerInfo().GetName().strip() == anion_oxygen:

                print(atom.GetImplicitValence(), atom.GetExplicitValence())
                print("found anionic oxygen")

                anion_position = mol.GetConformer().GetAtomPosition(atom.GetIdx())

                anion_pos_coords = [anion_position.x, anion_position.y, anion_position.z]
                
                
                matches = mol.GetSubstructMatches(bare_nitrogen)
                mol = Chem.RWMol(mol)
                
                match_idx = 0
                for match in matches:
                    match_idx = match[0]
                print(match_idx)
                print(atom.GetMonomerInfo().GetName())
                # print('here2')
                # # monomer_info = AtomMonomerInfo().SetName("Na1")
                # print('here3')
                # # sodium_atom.SetMonomerInfo(monomer_info)
                # print('here4')

                # new_atom = mol.AddAtom(sodium_atom)
                # print('here5')
                mol.AddBond(match_idx, atom.GetIdx(), order=Chem.BondType.IONIC)
                # print('here6')
                new_atom_position =2/np.linalg.norm([0,0,-1])*np.array(anion_pos_coords)
                # print('here7')
                mol.GetConformer().SetAtomPosition(match_idx, Point3D(new_atom_position[0],new_atom_position[1],new_atom_position[2]))
                # print('added sodium')
                atom.SetFormalCharge(-1)
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)
                atom.UpdatePropertyCache()
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        mol.RemoveAtom(neighbor.GetIdx())
                        return mol
                # print(atom.GetExplicitValence())

        return mol
    center = [0,15,1]
    axis_direction = normalize([0,0,1])

    angle = 0
    
    mol = Chem.MolFromPDBFile(input_struct, removeHs=True)
    
    core_mol = Chem.MolFromSmarts("C-S-C-C-1O-C-C-C-C-1")
    core_mol = Chem.RemoveAllHs(core_mol)
    print(Chem.MolToSmiles(core_mol))
    subunit_struct = Chem.MolFromSmiles("OC1OC(CSCCC(=O)[O-])CC(O)C1O")
    subunit_struct = Chem.RemoveAllHs(subunit_struct)
    print(Chem.MolToSmiles(subunit_struct))
    mol = AllChem.AssignBondOrdersFromTemplate(subunit_struct, mol)
    mol = Chem.RWMol(mol)
    print(len(mol.GetAtoms()))
    print('laksdfj')
   
    # return
    centroid = rdMolTransforms.ComputeCentroid(mol.GetConformer())
    center = [centroid.x, centroid.y, centroid.z]
    edit_mol = Chem.RWMol(mol)
    atom_idx_matches_start =[]
    atom_idx_matches_end =[]
    
    prev_mol = edit_mol
    
    bond_start = "O1"
    bond_end = "C7"

    plane_atoms = ["O2","C9","C7"]
    plane_atom_coords = []
    plane_norm_vec = []
    ring_atoms = ["O2","C9","C7","C1","C2","C8"]
    ring_atom_coords = []
    ring_center = []
    axis_direction_atoms = ["C2","C9"]
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
    axis_direction = normalize(np.subtract(axis_direction, axis_offset*np.array(center)))
    center = np.add(np.array(ring_center), (-center_offset/np.linalg.norm(plane_norm_vec)*np.array(plane_norm_vec)))
    atom_count = 1
    for i in range(8):
        angle += (2*np.pi/8)
        
        mol = Chem.MolFromPDBFile(input_struct, removeHs=True)
        mol = AllChem.AssignBondOrdersFromTemplate(subunit_struct, mol)
        Chem.AssignStereochemistryFrom3D(mol)
        mol = add_sodium(mol)
        print('super done')
        # Chem.MolToMolFile(mol, output_path)
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            position = conf.GetAtomPosition(atom.GetIdx())
            new_position = rot_point_from_formula(center[0],center[1],center[2],axis_direction[0],axis_direction[1],axis_direction[2],position.x, position.y, position.z, angle)

            conf.SetAtomPosition(atom.GetIdx(), Point3D(new_position[0], new_position[1], new_position[2]))
        mols.append(mol)

    combined_mols = functools.reduce(Chem.CombineMols, mols)
    
    for atom in combined_mols.GetAtoms():
        if not atom.GetMonomerInfo():
            continue
        if atom.GetMonomerInfo().GetName().strip() == bond_start:
            atom_idx_matches_start.append(atom.GetIdx())
        elif atom.GetMonomerInfo().GetName().strip() == bond_end:
            atom_idx_matches_end.append(atom.GetIdx())

    # combined_mols = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromMolFile("/Users/adam/Downloads/inputs_for_molec_replac/sugammadex_protonated.mol"),combined_mols)
    combined_mols = Chem.RWMol(combined_mols)
    
    match_structs = combined_mols.GetSubstructMatches(core_mol)
   
    
    for i in range(len(atom_idx_matches_start) - 1):
        combined_mols.AddBond(atom_idx_matches_start[i], atom_idx_matches_end[i+1], Chem.rdchem.BondType.SINGLE)
        
    
    combined_mols.AddBond(atom_idx_matches_start[-1], atom_idx_matches_end[0], Chem.rdchem.BondType.SINGLE)
    for atom_idx in atom_idx_matches_end:
        combined_mols.GetAtomWithIdx(atom_idx).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    print('herererere')
    # combined_mols = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromMolFile("/Users/adam/Downloads/inputs_for_molec_replac/sugammadex_protonated.mol"),combined_mols)

    Chem.SanitizeMol(combined_mols)
    combined_mols = Chem.RemoveAllHs(combined_mols)
    combined_mols = Chem.AddHs(combined_mols, addCoords=True, addResidueInfo=True)
    print(combined_mols)
    
    # combined_mols = add_sodium(combined_mols)
    conf = combined_mols.GetConformer()
    # ff = rdForceFieldHelpers.UFFGetMoleculeForceField(combined_mols)
    print('herererere')
    Chem.SanitizeMol(combined_mols)
    mp =  rdForceFieldHelpers.MMFFGetMoleculeProperties(combined_mols, mmffVerbosity=0)
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(combined_mols, mp)
    # for match_indices in match_structs:
    #     print(match_indices)
    #     for match_idx in match_indices:
    #         ff.MMFFAddPositionConstraint(match_idx, 0.5, 1.e4)
    #         # continue
    for atom in combined_mols.GetAtoms():
        # if atom.GetAtomicNum() in [1]:
        #     continue
        # if atom.GetAtomicNum() in [7]:
        #     ff.MMFFAddPositionConstraint(atom.GetIdx(), , 1.e4)
        ff.MMFFAddPositionConstraint(atom.GetIdx(), constraint, 1.e4)
    # rdForceFieldHelpers.MMFFOptimizeMolecule(combined_mols, maxIters=1000000)
    
    Chem.SanitizeMol(combined_mols)
    # ff.Minimize(maxIts=10000)
    # combined_mols = Chem.ReplaceSubstructs(combined_mols, Chem.MolFromSmiles("[NH4+]"), Chem.MolFromSmiles("[Na+]"), replaceAll=True)
    combined_mols = Chem.RWMol(combined_mols)
    n_query = rdqueries.AtomNumEqualsQueryAtom(7)
    are_N_present = len(combined_mols.GetAtomsMatchingQuery(n_query))
    
    count = 0
    e_brake = False
    while count < 32 and not e_brake:
        for atom in combined_mols.GetAtoms():
            has_broken = False
            if atom.GetAtomicNum() == 7:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        combined_mols.RemoveAtom(neighbor.GetIdx())
                        has_broken = True
                        count += 1
                        print(count)
                        break
            if has_broken:
                print('breaking')
                break
            # print('exiting while loop')
            # e_brake = True
    for atom in combined_mols.GetAtoms():
        if atom.GetAtomicNum() == 7:
            print('found dummy atom, replacing')
            atom.SetAtomicNum(11)
            atom.SetFormalCharge(1)
            # atom.SetNumExplicitHs(0)
            # atom.SetNoImplicit(True)
            atom.UpdatePropertyCache()
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    combined_mols.RemoveAtom(neighbor.GetIdx())
    while (len(combined_mols.GetAtomsMatchingQuery(n_query))) > 0:
        print('in loop')
        for atom in combined_mols.GetAtoms():
            if atom.GetAtomicNum() == 7:
                combined_mols.RemoveAtom(atom.GetIdx())
                break

    for atom in combined_mols.GetAtoms():
            atom_mon_inf = AtomMonomerInfo()
            atom_mon_inf.SetName(f"{atom.GetSymbol}{atom_count}")
            atom.SetMonomerInfo(atom_mon_inf)
            atom_count += 1
    # combined_mols = Chem.RemoveHs(combined_mols, updateExplicitCount=True)
    # Chem.SanitizeMol(combined_mols)
    # combined_mols = Chem.AddHs(combined_mols, addCoords=True,addResidueInfo=True)
    print('succesfully replaced', combined_mols)
    # output_path = output_path.split(".pdb")[0] + str(axis_offset) + ".pdb"
    Chem.MolToPDBFile(combined_mols, output_path)
    []

    
if __name__ == "__main__":
    input_dir = "/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_2"
    output_dir = "/Users/adam/Downloads/inputs_for_molec_replac/SUG_CONFORGE_TRIAL_15"
    count = 0
    sorted_dir = []
    for f in os.listdir(input_dir):
        if ".pdb" not in f or "272" not in f:
            continue
        # for i in range(500):
        sorted_dir.append(f)
    sorted_dir = list(sorted(sorted_dir, key = lambda x : int(x.split("_")[-1].split(".pdb")[0])))
    # print(sorted_dir[104])
    # raise Exception('test')
    for f in sorted_dir:
        print(f)
        new_name = f.split(".pdb")[0]
        try:
            symmetrize(f"{input_dir}/{f}",f"{output_dir}/{new_name}_symmetrized_{count}.pdb", constraint=0.0, center_offset=4.9, axis_offset=0)
        except Exception as e:
            print(e)
            continue
        count += 1  