# import gemmi
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem, rdForceFieldHelpers, rdchem, rdmolops, Draw
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
input_path = "/Users/adam/Downloads/inputs_for_molec_replac/simeprevir.pdb"


def enumerateTorsions(mol):
   torsionSmarts = '[!$(*#*)&!D1]~[!$(*#*)&!D1]'
   torsionQuery = Chem.MolFromSmarts(torsionSmarts)
   matches = mol.GetSubstructMatches(torsionQuery)
   torsionList = []
   for match in matches:
     idx2 = match[0]
     idx3 = match[1]
     bond = mol.GetBondBetweenAtoms(idx2, idx3)
     jAtom = mol.GetAtomWithIdx(idx2)
     kAtom = mol.GetAtomWithIdx(idx3)
     if (((jAtom.GetHybridization() != Chem.HybridizationType.SP2)
       and (jAtom.GetHybridization() != Chem.HybridizationType.SP3))
       or ((kAtom.GetHybridization() != Chem.HybridizationType.SP2)
       and (kAtom.GetHybridization() != Chem.HybridizationType.SP3))):
       continue
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
         torsionList.append((idx1, idx2, idx3, idx4))
   return torsionList
mol = Chem.MolFromPDBFile(input_path)
dihed_angles = enumerateTorsions(mol)
allowed_dihed_angles = []

for angle in dihed_angles:
    try:
        orig_angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(),angle[0],angle[1],angle[2],angle[3])
        rdMolTransforms.SetDihedralDeg(mol.GetConformer(),angle[0],angle[1],angle[2],angle[3], orig_angle)
    except:
       continue
    allowed_dihed_angles.append(angle)
allowed_dihed_angle_names = []
for allowed_angle in allowed_dihed_angles:
   print(mol.GetAtomWithIdx(allowed_angle[0]).GetMonomerInfo().GetName())
print(allowed_dihed_angles)

