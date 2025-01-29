from Bio.PDB import PDBParser, Superimposer, Selection
from rdkit.Chem import rdMolAlign, AllChem
from rdkit import Chem
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


if __name__ == "__main__":
    print(calc_rmsd("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha_reordered_2.pdb","/Users/adam/Downloads/outputs_from_molec_replac/PAR_CUSTOM_CONF_TRIAL_2_15A/ROUND_10/paritaprevir_torsion_angle_perturb_32_out.1_refine_001.pdb"))