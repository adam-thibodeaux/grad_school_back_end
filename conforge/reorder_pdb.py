import numpy as np
import rdkit.Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdMolAlign

def reorder_atoms(mol_pdb_fname, template_pdb_fname, output_pdb_fname):
    from rdkit.Chem import rdmolfiles
    bond_template = rdkit.Chem.MolFromMolFile("/Users/adam/Downloads/paritaprevir_correct_bonds.mol",removeHs=False)
    
    
    mol_to_transform = rdkit.Chem.rdmolfiles.MolFromPDBFile(mol_pdb_fname, removeHs=False)
    AllChem.AssignBondOrdersFromTemplate(bond_template,mol_to_transform)
    transform_order = list(rdmolfiles.CanonicalRankAtoms(mol_to_transform))

    mol_template = rdkit.Chem.rdmolfiles.MolFromPDBFile(template_pdb_fname, removeHs=False)
    AllChem.AssignBondOrdersFromTemplate(bond_template,mol_template)
    template_order = list(rdmolfiles.CanonicalRankAtoms(mol_template))

    if len(template_order) != len(transform_order):
        raise RuntimeError('Number of atoms differs between template and molecule to transform.')

    i_transform_order = [int(i) for i in np.argsort(transform_order)]
    i_template_order =  [int(i) for i in np.argsort(template_order)]

    N_atoms = len(template_order)

    pos_to_transform = mol_to_transform.GetConformers()[0].GetPositions()
    lines = [None]*N_atoms
    for _, (otr, ote) in enumerate(zip(i_transform_order, i_template_order)):
        print(mol_to_transform.GetAtoms()[otr].GetPDBResidueInfo().GetName(),
             mol_template.GetAtoms()[ote].GetPDBResidueInfo().GetName())
        pdb_entry_template = mol_template.GetAtoms()[ote].GetPDBResidueInfo()
        lines[ote] = '{ATOM:<6}{serial_number:>5} {atom_name:<4}{alt_loc_indicator:<1}{res_name:<3} {chain_id:<1}{res_seq_number:>4}{insert_code:<1}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}'.format(
                    ATOM='HETATM',
                    serial_number=pdb_entry_template.GetSerialNumber(),
                    atom_name=pdb_entry_template.GetName(),
                    alt_loc_indicator=' ',
                    res_name=pdb_entry_template.GetResidueName(),
                    chain_id=pdb_entry_template.GetChainId(),
                    res_seq_number=pdb_entry_template.GetResidueNumber(),
                    insert_code=' ',
                    x= pos_to_transform[otr, 0],
                    y= pos_to_transform[otr, 1],
                    z= pos_to_transform[otr, 2],
                    occupancy=1.0,
                    temp_factor=0.0
                )

    with open(output_pdb_fname, 'w') as f:
        for line in lines:
            if line is not None:
                print(line, file=f)


def validate(template_pdb_fname, output_pdb_fname):
    import rdkit.Chem.rdPartialCharges
    mol_template = rdkit.Chem.rdmolfiles.MolFromPDBFile(template_pdb_fname, removeHs=False)
    finished_molecule = rdkit.Chem.rdmolfiles.MolFromPDBFile(output_pdb_fname, removeHs=False)
    N_atoms = finished_molecule.GetNumAtoms()
    for i in range(N_atoms):
        assert finished_molecule.GetAtoms()[i].GetAtomicNum() == mol_template.GetAtoms()[i].GetAtomicNum()
        assert finished_molecule.GetAtoms()[i].GetDegree() == mol_template.GetAtoms()[i].GetDegree()
        assert finished_molecule.GetAtoms()[i].GetTotalDegree() == mol_template.GetAtoms()[i].GetTotalDegree(), i
        assert finished_molecule.GetAtoms()[i].GetHybridization() == mol_template.GetAtoms()[i].GetHybridization(), i
        assert finished_molecule.GetAtoms()[i].GetFormalCharge() == mol_template.GetAtoms()[i].GetFormalCharge(), i
        assert finished_molecule.GetAtoms()[i].GetTotalValence() == mol_template.GetAtoms()[i].GetTotalValence(), i    
    rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(finished_molecule, throwOnParamFailure=True)
    charges_finished = np.array([float(finished_molecule.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(N_atoms)])
    rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol_template, throwOnParamFailure=True)
    charges_template = np.array([float(mol_template.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(N_atoms)])
    assert np.allclose(charges_finished, charges_template)
    print('Validation OK')

def reorder_atoms_v2(input_path, reference_path, output_path):
    input_pdb = rdkit.Chem.MolFromPDBFile(input_path,removeHs=False)
    reference_pdb = rdkit.Chem.MolFromPDBFile(reference_path,removeHs=False)
    print(len(input_pdb.GetAtoms()))
    print(len(reference_pdb.GetAtoms()))
    # input_pdb_atom_dict = {}
    # reference_pdb_atom_dict = {}
    # mapping = []
    # for atom in input_pdb.GetAtoms():
    #     print(atom.GetMonomerInfo())
    #     input_pdb_atom_dict[atom.GetMonomerInfo().GetName().strip()] = atom.GetIdx()
    # for atom in reference_pdb.GetAtoms():
    #     mapping.append( (atom.GetIdx(),input_pdb_atom_dict[atom.GetMonomerInfo().GetName().strip()]))
    # rdkit.Chem.rdMolAlign.AlignMol(input_pdb,reference_pdb)
    atoms = input_pdb.GetSubstructMatch(reference_pdb)
    print(atoms)
    new_pdb = rdkit.Chem.RenumberAtoms(input_pdb, atoms)
    rdkit.Chem.MolToPDBFile(new_pdb, output_path)

if __name__ == '__main__':
    reorder_atoms(
       "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta.pdb",
       "/Users/adam/Downloads/outputs_from_molec_replac/PAR_BETA_CUSTOM_CONF_TRIAL_1/ROUND_10/paritaprevir_torsion_angle_perturb_19_out.1_refine_001.pdb",
        "/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_beta_reordered.pdb",
                  )
    # validate("/Users/adam/Downloads/inputs_for_molec_replac/PAR_CONFORGE_TRIAL_1/paritaprevir_alpha_conforge_78.pdb","/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha_reordered_2.pdb")