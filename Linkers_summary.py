# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:54:38 2024

@author: Coretshz
"""

import shutil
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

cycle_number=1
os.makedirs('Linker_summary/' + 'Cycle ' + str(cycle_number))

for i in os.listdir('optimal_linker'):
    shutil.copy('optimal_linker/' + i, 'Linker_summary/Cycle ' + str(cycle_number) + '/' + i)
    mol = Chem.MolFromMolFile('Linker_summary/Cycle ' + str(cycle_number) + '/' + i)
    if mol is None: continue
    Chem.SanitizeMol(mol)
    mol_r_Hs = Chem.AddHs(mol) # Adding H atoms
    AllChem.EmbedMolecule(mol_r_Hs, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol_r_Hs)
    Chem.MolToMolFile(mol_r_Hs, filename = 'Linker_summary/Cycle ' + str(cycle_number) + '/' + i)
    AllChem.Compute2DCoords(mol_r_Hs)
    Draw.MolToFile(mol_r_Hs, molsPerRow = 4, subImgSize = (200, 200), filename = 'Linker_summary/Cycle ' + str(cycle_number) + '/' + i.split('.')[0] + '.png')