# -*- coding: utf-8 -*-

def SL(optimal_linker, sdf_name, best_edges_png, best_edges_mol):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
    from sklearn.metrics import jaccard_score
    import numpy as np
    
    suppl = Chem.SDMolSupplier(sdf_name)
    mols1 = [Y for Y in suppl if Y is not None]
    len(mols1)
    maccs_fps_500000 = [AllChem.GetMACCSKeysFingerprint(mol1) for mol1 in mols1]
    fps_si = []
    best_si = []
    best_si_index = []
    for Q in maccs_fps_500000:
        S = jaccard_score(optimal_linker, np.array(Q))
        fps_si.append(S)
        if S > 0.8:
            best_si.append(S)
            best_si_index.append(maccs_fps_500000.index(Q))

    best_edges = []
    for R in best_si_index:
        mol_r = mols1[R]
        best_edges.append(mol_r)
        mol_r_Hs = Chem.AddHs(mol_r) # Outputing CIF file
        AllChem.EmbedMolecule(mol_r_Hs)
        AllChem.UFFOptimizeMolecule(mol_r_Hs)
        Chem.MolToMolFile(mol_r_Hs, filename = best_edges_mol + '/' + str(R) + '.mol')
        AllChem.Compute2DCoords(mol_r_Hs)
        Draw.MolToFile(mol_r_Hs, molsPerRow = 4, subImgSize = (200, 200), filename = best_edges_png + '/' + str(R) + '.png')