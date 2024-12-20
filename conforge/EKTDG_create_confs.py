from rdkit import Chem
from rdkit.Chem import rdDistGeom, AllChem, TorsionFingerprints
from rdkit.Geometry.rdGeometry import Point3D
import subprocess
import cpeptools
import cpeptools.cpeptools
import os
def create_conformers(is_restricted=False):
    # need this to ensure that double bonds are maintainied and right geometry
    mol_template = Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds.mol", removeHs=False)
    # Chem.MolToPDBFile(mol, "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha_no_H.pdb")
    # return


    mol = Chem.MolFromPDBFile("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb", removeHs=False, sanitize=True, proximityBonding=False)
    # print(mol_template)
    mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
    if is_restricted:
        restricted_frag = Chem.MolFromPDBFile("/Users/adam/Downloads/inputs_for_molec_replac/PAR_FRAG_EKTDG_TRIAL_1/ROUND_1_2/paritaprevir_core_only_ektdg_79.pdb", sanitize=False, removeHs=True)
    # print(restricted_frag)
        restricted_conf = restricted_frag.GetConformer()
        restricted_atoms = {}
        for atom in restricted_frag.GetAtoms():
            # print(atom.GetMonomerInfo().GetName())
            coords = restricted_conf.GetAtomPosition(atom.GetIdx())
            atom_coord = Point3D()
            atom_coord.x = coords.x
            atom_coord.y = coords.y
            atom_coord.z = coords.z
            restricted_atoms[atom.GetIdx()] = atom_coord
    
    # return
    # print(restricted_atoms)
    ps = rdDistGeom.ETKDGv3()
    ps.pruneRmsThresh = 0.1
    ps.numThreads = 0
    # ps.useBasicKnowledge = True
    # ps.useRandomCoords = True
    # ps.optimizerForceTol = 0.0005
    ps.randomSeed = 210185

    # ps.useMacrocycle14config = False
    # ps.useMacrocycleTorsions = False
    # ps.useExpTorsionAnglePrefs = False

    # bmat = cpeptools.cpeptools.bound_matrix_from_ellipse(mol, angle = 0, eccentricity = 0.99)
    # ps.SetBoundsMat(bmat)
    if is_restricted: 
        ps.SetCoordMap(restricted_atoms)
    # ps.ignoreSmoothingFailures = True
    # ps.ETversion = 3
    # ps.coordMap = restricted_atoms
    cids = rdDistGeom.EmbedMultipleConfs(mol,300,ps)
    count = 0

    # return
    print(len(cids))
    tfd_vals = TorsionFingerprints.GetTFDBetweenConformers(mol,[0],[i for i in range(1,len(cids))], useWeights=False)
    tfd_vals_conf_map = {}
    for i in range(len(tfd_vals)):
        tfd_val = tfd_vals[i]
        print(tfd_val)
        tfd_vals_conf_map[i] = tfd_val
    for cid in cids:
        # mol = Chem.AddHs(mol)
        Chem.MolToPDBFile(mol, f"/Users/adam/Downloads/inputs_for_molec_replac/PAR_BETA_EKTDG/paritaprevir_beta_ektdg_{count}.pdb", cid)
        count += 1
    with open("/Users/adam/Downloads/inputs_for_molec_replac/PAR_EKTDG_TRIAL_1_with_H/tfd_out.txt", "w") as file:
        file.write(str(sorted(tfd_vals_conf_map.items(), key=lambda x : x[1])))
        file.close()

            

def correlate_csd_pdb():
    # input_path = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb"
    # output_path = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha_solved_EDKTG_paired.pdb"
    # mol = Chem.MolFromPDBFile(input_path, removeHs=False)
    string = """
HETATM    1  C1  UNK     1       1.914  -2.249  45.611  1.00  1.11           C
HETATM    2  C2  UNK     1       1.624  -1.044  46.428  1.00  1.42           C
HETATM    3  C3  UNK     1       1.099   0.012  45.413  1.00  0.39           C
HETATM    4  C4  UNK     1       0.407  -0.841  44.316  1.00  0.79           C
HETATM    5  C5  UNK     1       0.652  -0.303  42.884  1.00  0.32           C
HETATM    6  C6  UNK     1      -0.204   0.606  40.771  1.00  0.79           C
HETATM    7  C7  UNK     1      -1.379   0.322  39.898  1.00  0.32           C
HETATM    8  C8  UNK     1      -0.214  -0.638  39.827  1.00  0.87           C
HETATM    9  C9  UNK     1       0.708  -0.816  38.598  1.00  1.42           C
HETATM   10  C10 UNK     1       1.227  -2.028  38.192  1.00  0.95           C
HETATM   11  C11 UNK     1       0.931  -3.337  38.771  1.00  1.11           C
HETATM   12  C12 UNK     1       2.138  -4.007  39.446  1.00  0.63           C
HETATM   13  C13 UNK     1       1.812  -5.399  40.050  1.00  0.39           C
HETATM   14  C14 UNK     1       0.972  -5.331  41.355  1.00  1.03           C
HETATM   15  C15 UNK     1       1.863  -4.955  42.569  1.00  1.50           C
HETATM   16  C16 UNK     1       1.110  -4.572  43.838  1.00  0.39           C
HETATM   17  C18 UNK     1      -0.723  -5.549  45.321  1.00  0.63           C
HETATM   18  C17 UNK     1       0.382  -3.149  43.706  1.00  0.71           C
HETATM   19  C19 UNK     1      -1.680  -6.680  45.468  1.00  1.26           C
HETATM   20  C20 UNK     1      -2.632  -6.706  46.494  1.00  1.58           C
HETATM   21  C22 UNK     1      -4.408  -9.897  46.037  1.00  2.37           C
HETATM   22  C21 UNK     1      -3.502  -8.731  45.798  1.00  1.26           C
HETATM   23  C23 UNK     1      -2.504  -8.709  44.722  1.00  1.58           C
HETATM   24  C24 UNK     1       0.596   1.761  40.462  1.00  0.79           C
HETATM   25  C25 UNK     1       1.084   5.117  40.685  1.00  1.26           C
HETATM   26  C26 UNK     1       1.496   6.583  40.954  1.00  1.82           C
HETATM   27  C27 UNK     1       0.188   5.980  41.538  1.00  1.42           C
HETATM   28  C28 UNK     1       0.295  -0.620  48.475  1.00  0.79           C
HETATM   29  C29 UNK     1      -0.728  -1.186  49.399  1.00  1.66           C
HETATM   30  C30 UNK     1      -1.446  -2.334  49.018  1.00  1.03           C
HETATM   31  C31 UNK     1      -2.423  -2.763  49.927  1.00  1.42           C
HETATM   32  C32 UNK     1      -2.672  -2.146  51.130  1.00  1.11           C
HETATM   33  C33 UNK     1      -1.853  -0.988  51.491  1.00  0.39           C
HETATM   34  C34 UNK     1      -0.911  -0.517  50.592  1.00  1.03           C
HETATM   35  C35 UNK     1      -0.056   0.662  50.958  1.00  1.11           C
HETATM   36  C36 UNK     1      -0.143   1.474  52.090  1.00  1.50           C
HETATM   37  C37 UNK     1       0.723   2.468  52.364  1.00  1.11           C
HETATM   38  C38 UNK     1       1.731   2.844  51.394  1.00  1.11           C
HETATM   39  C39 UNK     1       1.776   2.118  50.232  1.00  1.34           C
HETATM   40  C40 UNK     1       0.947   1.046  49.942  1.00  0.63           C
HETATM   84  N1  UNK     1       0.937  -2.159  44.499  1.00  0.16           N
HETATM   85  N2  UNK     1      -0.372   0.223  42.213  1.00  0.55           N
HETATM   86  N3  UNK     1       1.074   2.569  41.563  1.00  1.11           N
HETATM   87  N4  UNK     1       0.081  -5.635  44.158  1.00  0.87           N
HETATM   88  N5  UNK     1      -1.685  -7.694  44.600  1.00  1.18           N
HETATM   89  N6  UNK     1      -3.543  -7.665  46.692  1.00  1.11           N
HETATM   90  N7  UNK     1       1.079   0.375  48.739  1.00  0.63           N
HETATM   91  O1  UNK     1      -0.529  -2.969  42.940  1.00  1.50           O
HETATM   92  O2  UNK     1       1.832  -0.367  42.452  1.00  1.03           O
HETATM   93  O3  UNK     1       0.916   2.156  39.339  1.00  0.71           O
HETATM   94  O4  UNK     1       3.232   3.518  40.720  1.00  1.42           O
HETATM   95  O5  UNK     1       2.209   4.262  42.868  1.00  2.05           O
HETATM   96  O6  UNK     1      -0.626  -4.619  46.129  1.00  1.50           O
HETATM   97  O7  UNK     1       0.443  -1.374  47.291  1.00  1.03           O
HETATM   98  S1  UNK     1       2.033   3.862  41.514  1.00  1.17           S
HETATM   41  H1  UNK     1       1.780  -3.059  46.128  1.00  1.34           H
HETATM   49  H2  UNK     1       2.400  -0.731  46.940  1.00  1.66           H
HETATM   56  H3  UNK     1       0.471   0.618  45.832  1.00 18.16           H
HETATM   70  H4  UNK     1      -0.560  -0.860  44.490  1.00  0.95           H
HETATM   81  H5  UNK     1      -1.521   0.906  39.137  1.00  0.39           H
HETATM   82  H6  UNK     1      -0.380  -1.478  40.306  1.00  1.03           H
HETATM   83  H7  UNK     1       0.925  -0.060  38.102  1.00  1.66           H
HETATM   42  H8  UNK     1       1.819  -2.010  37.474  1.00  1.11           H
HETATM   43  H9  UNK     1       0.601  -3.922  38.073  1.00  1.34           H
HETATM   44  H10 UNK     1       2.465  -3.427  40.148  1.00  0.79           H
HETATM   45  H11 UNK     1       2.644  -5.861  40.237  1.00  0.47           H
HETATM   46  H12 UNK     1       0.268  -4.671  41.255  1.00  1.26           H
HETATM   47  H13 UNK     1       2.443  -5.705  42.767  1.00  1.82           H
HETATM   48  H14 UNK     1       1.751  -4.528  44.579  1.00  0.47           H
HETATM   50  H15 UNK     1      -2.630  -5.994  47.094  1.00  1.89           H
HETATM   51  H16 UNK     1      -4.590 -10.336  45.202  1.00  3.63           H
HETATM   52  H17 UNK     1      -2.452  -9.419  44.124  1.00  1.89           H
HETATM   53  H18 UNK     1       0.768   4.905  39.784  1.00  1.50           H
HETATM   54  H19 UNK     1       1.417   7.216  40.221  1.00  2.21           H
HETATM   55  H20 UNK     1       0.160   5.801  42.491  1.00  1.66           H
HETATM   57  H21 UNK     1      -1.282  -2.779  48.219  1.00  1.18           H
HETATM   58  H22 UNK     1      -2.936  -3.506  49.699  1.00  1.74           H
HETATM   59  H23 UNK     1      -3.341  -2.456  51.700  1.00  1.34           H
HETATM   60  H24 UNK     1      -1.965  -0.574  52.317  1.00  0.47           H
HETATM   61  H25 UNK     1       2.409   2.361  49.595  1.00  1.58           H
HETATM   62  H26 UNK     1       2.317   3.546  51.555  1.00  1.34           H
HETATM   63  H27 UNK     1       0.670   2.914  53.179  1.00  1.26           H
HETATM   64  H28 UNK     1      -0.837   1.317  52.688  1.00  1.74           H
HETATM   65  H29 UNK     1       2.825  -2.233  45.279  1.00  1.34           H
HETATM   66  H30 UNK     1      -1.124   0.344  42.609  1.00  0.71           H
HETATM   67  H31 UNK     1       1.832   0.528  45.037  1.00  9.47           H
HETATM   68  H32 UNK     1       0.575   2.185  42.447  1.00  1.66           H
HETATM   69  H33 UNK     1      -0.193  -6.260  43.554  1.00  1.26           H
HETATM   71  H34 UNK     1      -2.189   0.010  40.326  1.00  0.39           H
HETATM   72  H35 UNK     1       0.227  -3.235  39.430  1.00  1.34           H
HETATM   73  H36 UNK     1       2.844  -4.107  38.789  1.00  0.79           H
HETATM   74  H37 UNK     1       1.326  -5.922  39.393  1.00  0.47           H
HETATM   75  H38 UNK     1       0.555  -6.192  41.519  1.00  1.26           H
HETATM   76  H39 UNK     1       2.429  -4.209  42.311  1.00  1.82           H
HETATM   77  H40 UNK     1      -5.231  -9.588  46.422  1.00  3.63           H
HETATM   78  H41 UNK     1      -3.985 -10.514  46.637  1.00  3.63           H
HETATM   79  H42 UNK     1       2.233   6.745  41.560  1.00  2.21           H
HETATM   80  H43 UNK     1      -0.655   6.270  41.154  1.00  1.66           H
"""
    ordered_string = string.split("\n")
    # ordered_string = sorted(split_string, key= lambda x: x[13:15])
    print(ordered_string)
    for i in range(len(ordered_string)):
        line = ordered_string[i]
        if i < 10:
            ordered_string[i] = f"{line[:9]} {i}{line[11:]}".replace("UNK", "UNL").strip()
        else:
            ordered_string[i] = f"{line[:9]}{i}{line[11:]}".replace("UNK", "UNL").strip()
    print("\n".join(ordered_string))

def add_hydrogens(input_path, output_path):
    
    # babel_str = f"obabel {input_path} -opdb {output_path} -h"
    # print(babel_str)
    # subprocess.run([babel_str])
    mol = Chem.MolFromPDBFile(input_path, sanitize=False, removeHs=True)
    mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToPDBFile(mol, output_path)
def renumber_atoms(input_path, reference_path, output_path):
    input_struct = Chem.MolFromPDBFile(input_path, removeHs=False)
    reference_struct = Chem.MolFromPDBFile(reference_path, removeHs=False)
    match = input_struct.GetSubstructMatch(reference_struct)
    Chem.MolToPDBFile(Chem.RenumberAtoms(input_struct, match), output_path)
    print(match)
if __name__ == "__main__":
    create_conformers(is_restricted=False)
    # correlate_csd_pdb()
    # add_hydrogens("/Users/adam/Downloads/inputs_for_molec_replac/PAR_EKTDG_TRIAL_1/paritaprevir_ektdg_5345.pdb", "/Users/adam/Downloads/inputs_for_molec_replac/PAR_EKTDG_TRIAL_1/_with_H_paritaprevir_ektdg_5345.pdb")
    # for f in os.listdir("/Users/adam/Documents/code/grad_school/back_end/conforge/para_conformers_1000_conf_120_e_0_01_rmsd/"):
    #     input_path = f"/Users/adam/Documents/code/grad_school/back_end/conforge/para_conformers_1000_conf_120_e_0_01_rmsd/{f}"
    #     reference_path = "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.pdb"
    #     output_path = f"/Users/adam/Downloads/inputs_for_molec_replac/PAR_CONFORGE_120_E_0_01_RMSD/{f}"
    #     # add_hydrogens(input_path=input_path, output_path=output_path)
    #     try:
    #         renumber_atoms(input_path=input_path, reference_path=reference_path, output_path=output_path)
    #     except Exception as e:
    #         print(e)
    #         pass
