import rdkit.Chem as Chem
from rdkit.Chem import rdForceFieldHelpers, AllChem
import numpy as np
import matplotlib.pyplot as plt
from rdkit.Geometry import Point3D
import vg
import gemmi
input_struct = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta.pdb"
# atom_names = ["C11","C12","C13", "C14"]
atom_names = ["C14","C13","C12","C11"]
adjacent_atom_names = ["H12A","H12B"]
output_path = "/Users/adam/Downloads/dihedral_change/"

def calc_new_fit(pdb_path, mtz_path, count, orig_pdb):
    try:
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
            intensity = ccp4.grid.interpolate_value(gemmi.Position(position.x,position.y,position.z))
            orig_position = orig_mol.GetConformer(0).GetAtomPosition(i)
            if abs(orig_position.x - position.x) < 0.01 and abs(orig_position.y - position.y) < 0.01 and abs(orig_position.z - position.z) < 0.01:
                i+= 1
                continue
            if intensity < 0:
                intensities.append(abs(intensity))
            i+= 1
        # position = this_mol.GetConformer(cid).GetAtomPosition(atom_idx_input)
        #                 # print(position.x - base_position.x, position.y-base_position.y, position.z-base_position.z)
        # intensity = ccp4.grid.interpolate_value(gemmi.Position(position.x,position.y,position.z))
        # print(sum(intensities))
        return [count[0], sum(intensities)]
    except Exception as e:
        print(e)
        return 

def circularpoints_v3(center, radius, normal, num_points=12, offset=0):
# normal == direction Vector (V)

    v1 =v2 = np.array([0. , 0., 0.]).astype(np.float64) #z
    if True:
        if abs(normal[0]) > abs(normal[2]):
            b2 = np.array([-normal[1], normal[0], 0.0])
        else:
            b2 = np.array([0.0, -normal[2], normal[1]])
        b2 /= np.linalg.norm(b2)
        b1 = np.cross(b2, normal)
        b1 /= np.linalg.norm(b1)
        v1 = b2
        v2 = b1
    circle = []
    angles = np.linspace(0 + offset, 2 * np.pi + offset, num=num_points)
    for angle in angles:
        v1_Out = np.cos(angle) * v1 
        v2_Out = np.sin(angle) * v2
        point = center + (radius * (v1_Out + v2_Out ))
        circle.append(point)
    # print(circle)
    return circle
def main(input_struct=input_struct,atom_names=atom_names,adjacent_atom_names = adjacent_atom_names, output_path=output_path, num_points = 12, count = [0,[]], mtz_path= ""):
    mol = Chem.MolFromPDBFile(input_struct, removeHs=False, sanitize=False)
    # template_mol = Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds.mol", removeHs=False)
    # mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    print(atom_names)
    print(adjacent_atom_names)
    atom_idx_dict = {}
    phaser_mtz = input_struct.split("/")
    name = phaser_mtz.pop()
    phaser_mtz.append(name.split("_out")[0]+"_out.1.mtz")
    phaser_mtz = "/".join(phaser_mtz)
    for i, atom in enumerate(mol.GetAtoms()):
        atom_name = atom.GetMonomerInfo().GetName().strip()
        atom_idx_dict[atom_name] = atom.GetIdx()
    for i in range(0,1):
            conformer = mol.GetConformer()
            super_left_idx = atom_idx_dict[atom_names[0]]
            left_idx = atom_idx_dict[atom_names[1]]
            right_idx = atom_idx_dict[atom_names[2]]
            super_right_idx = atom_idx_dict[atom_names[3]]
            coords_super_left = conformer.GetAtomPosition(super_left_idx)
            coords_left = conformer.GetAtomPosition(left_idx) 
            coords_right = conformer.GetAtomPosition(right_idx) 
            coords_super_right = conformer.GetAtomPosition(super_right_idx) 
            np_coords_super_left = np.array([coords_super_left.x,coords_super_left.y,coords_super_left.z])
            np_coords_left = np.array([coords_left.x,coords_left.y,coords_left.z])
            np_coords_right = np.array([coords_right.x,coords_right.y,coords_right.z])
            np_coords_super_right = np.array([coords_super_right.x, coords_super_right.y,coords_super_right.z])
            center = (np_coords_left + np_coords_super_right)/2
            distance = np.sqrt(np.sum(np.square(center - np_coords_right)))
            adj_new_coords = []
            for adj_atom in adjacent_atom_names:
                adj_idx = atom_idx_dict[adj_atom]
                adj_coords = conformer.GetAtomPosition(adj_idx)
                np_adj_coords = np.array([adj_coords.x, adj_coords.y, adj_coords.z])
                adj_distance = np.sqrt(np.sum(np.square(center - np_adj_coords)))
                # offset_cos = np.dot(center - np_adj_coords, center - np_coords_left)/(np.linalg.norm(center - np_adj_coords) * np.linalg.norm(center - np_coords_left))
                # offset = np.arccos(np.clip(offset_cos, -1, 1))
                
                offset = vg.signed_angle(center - np_coords_right,center - np_adj_coords - np_coords_right,look=np_coords_right,units="rad")
                norm_offset = vg.signed_angle(center - np_coords_left, center, look=np_coords_left, units = "rad")

                adj_new_coords.append((adj_idx, circularpoints_v3(center, adj_distance, np_coords_left - np_coords_super_right, offset=offset, num_points=num_points)))
            # super_left_left = np_coords_super_left - np_coords_left
            # right_left = np_coords_right - np_coords_left
            ##find center (hard)
            #find radius, should be easier?
            #what is the normal then?
            coords_center = circularpoints_v3(center, distance, np_coords_left - np_coords_super_right, num_points=num_points)
            intensities = []
            for j in range(len(coords_center)):
                coord = coords_center[j]
                # print(len(coords_center))

                point = Point3D(coord[0],coord[1],coord[2])
                conformer.SetAtomPosition(right_idx,point)
                for adj_coord in adj_new_coords:
                    adj_point = Point3D(adj_coord[1][j][0],adj_coord[1][j][1],adj_coord[1][j][2])
                    conformer.SetAtomPosition(adj_coord[0],adj_point)
                count[0] += 1
                # Chem.SanitizeMol(mol)
                # print(len(mol.GetBonds()))
                # print(len(template_mol.))
                # if len(mol.GetBonds()) != len(template_mol.GetBonds()):
                #     print("invalid bonding found, skipping this one")
                #     continue
                Chem.MolToPDBFile(mol, f"{output_path}paritaprevir_perturb_custom_dihed_{count[0]}.pdb")
                
                new_intense = calc_new_fit(f"{output_path}paritaprevir_perturb_custom_dihed_{count[0]}.pdb", phaser_mtz, count, input_struct)

                if new_intense:
                    intensities.append(new_intense)

                #adjust bonded atoms to have same orientation as before
    intensities = list(filter(lambda x : x[1] != 0, intensities))
    sorted_intensities = sorted(intensities, key= lambda x : abs(x[1] - 0))
    count[1] += sorted_intensities
    count[1] = [x for x in count[1] if len(x) > 1]
    count[1] = sorted(count[1], key = lambda x: abs(x[1]-0))
    return count
            # cosine = np.dot(super_left_left,right_left)/(np.linalg.norm(super_left_left) * np.linalg.norm(right_left))
            # print(coords_left.x,coords_left.y,coords_left.z)
            # print(coords_right.x,coords_right.y,coords_right.z)
            # angle = np.arccos(cosine)
            # # print(np.degrees(angle))
            # print(angle)
            # # print(np.cos(angle) * np.sin(angle))
            # distance = np.linalg.norm(np_coords_right-np_coords_left)
            # print(distance)

    
if __name__ == "__main__":
    main()

#eq of cone at hkl, with angle 2theta
#(x - h)² / cos²θ + (y - k)² / cos²θ - (z - l)² / sin²θ = 0;

# import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.linalg import norm
# import pylab as plt

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection='3d')


# def truncated_cone(p0, p1, R0, R1, color):
#     """
#     Based on https://stackoverflow.com/a/39823124/190597 (astrokeat)
#     """
#     # vector in direction of axis
#     v = p1 - p0
#     # find magnitude of vector
#     mag = norm(v)
#     # unit vector in direction of axis
#     v = v / mag
#     # make some vector not in the same direction as v
#     not_v = np.array([1, 1, 0])
#     if (v == not_v).all():
#         not_v = np.array([0, 1, 0])
#     # make vector perpendicular to v
#     n1 = np.cross(v, not_v)
#     # print n1,'\t',norm(n1)
#     # normalize n1
#     n1 /= norm(n1)
#     # make unit vector perpendicular to v and n1
#     n2 = np.cross(v, n1)
#     # surface ranges over t from 0 to length of axis and 0 to 2*pi
#     n = 80
#     t = np.linspace(0, mag, n)
#     theta = np.linspace(0, 2 * np.pi, n)
#     # use meshgrid to make 2d arrays
#     t, theta = np.meshgrid(t, theta)
#     R = np.linspace(R0, R1, n)
#     # generate coordinates for surface
#     X, Y, Z = [p0[i] + v[i] * t + R *
#                np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
#     ax.plot_surface(X, Y, Z, color=color, linewidth=0, antialiased=False)


# A0 = np.array([1, 3, 2])
# A1 = np.array([8, 5, 9])
# ax.set_xlim(0, 10)
# ax.set_ylim(0, 10)
# ax.set_zlim(0, 10)
# truncated_cone(A0, A1, 1, 5, 'blue')
# plt.show()

# from scipy.optimize import fsolve

# def surface1(x, y, z):
#     return (7*np.sqrt(3)+12)*x**2+(7*np.sqrt(3)+6)* y**2+(7*np.sqrt(3)-4)*z**2-8*x*y-24*y*z-12*x*z


# def surface2(x, y, z):
#     return x**2 +y**2 +z**2 -1.57

# def equations(p):
#     x, y, z = p
#     return (surface1(x, y, z), surface2(x, y,z ))

# intersection_points = fsolve(equations, (0, 0, 0))

# print(intersection_points)

# import sympy as sp

# x, y, z = sp.symbols('x y z')

# # Define your surfaces
# surface1 = (7*np.sqrt(3)+12)*x**2+(7*np.sqrt(3)+6)* y**2+(7*np.sqrt(3)-4)*z**2-8*x*y-24*y*z-12*x*z
# surface2 = x**2 +y**2 +z**2 -1.57

# # Solve for the intersection
# intersection = sp.solve([surface1, surface2], [x, y, z])

# print(intersection)