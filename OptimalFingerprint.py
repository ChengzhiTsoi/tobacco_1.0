# -*- coding: utf-8 -*-

def OF(path):
    import os
    from rdkit import Chem
    import numpy as np
    from rdkit.Chem import AllChem
    
    olinker_dir = os.listdir(path)
    ol = []
    maccs_fps_ol2 = []
    maccs_fps_index = []
    
    for file_optimal_linker in olinker_dir: 
        MOL_optimal_linker = Chem.MolFromMolFile(path + '/' + file_optimal_linker)
        ol.append(MOL_optimal_linker)
        if MOL_optimal_linker is None: continue
        maccs_fps_ol1 = AllChem.GetMACCSKeysFingerprint(MOL_optimal_linker)
        maccs_fps_index.append(file_optimal_linker)
        maccs_fps_ol2.append(np.array(maccs_fps_ol1))
        
    sum_maccs_fps_ol1 = sum(np.array(maccs_fps_ol2))
    sum_maccs_fps_ol2 = np.argsort(-sum_maccs_fps_ol1)
    
    AA = []
    for I in range(30):
        I = 1
        AA.append(I)
    for J in range(137):
        J = 0
        AA.append(J)
    fps_ol = zip(sum_maccs_fps_ol2, AA)
    fps_ol = sorted(fps_ol)
    A_new,target_fp = zip(*fps_ol)
    optimal_linker = target_fp
    return optimal_linker