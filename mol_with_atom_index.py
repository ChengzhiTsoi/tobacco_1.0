# -*- coding: utf-8 -*-

def MWAI(i):
   from rdkit import Chem
   from rdkit.Chem import AllChem
   atoms = i.GetNumAtoms()
   for idx in range(atoms):
       i.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',
                                      str(i.GetAtomWithIdx(idx).GetIdx()))
   return i