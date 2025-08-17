# -*- coding: utf-8 -*-

def DXA(path_mol, best_edges_png, best_edges_cif):
    
    import os
    from rdkit import Chem
    import itertools
    from rdkit.Chem.Draw import rdMolDraw2D
    from mol_with_atom_index import MWAI
    
    path_list = os.listdir(path_mol)
    
    # Judging molecular symmetry
    symmetric_total = []
    not_symmetric_total = []
    central_atom_total = []
    up_atom_total = []
    down_atom_total = []
    symmetrical_molecule = []
    not_symmetrical_molecule = []
    ee = 0
    ff = 0
    for filename in path_list:
        if '.mol' not in filename:
            pass
        else:
            e = open(os.path.join(path_mol, filename), 'rb')
            mol2 = Chem.MolFromMolFile(os.path.join(path_mol, filename))
            atoms2 = mol2.GetAtoms()
            at = []
            atom_symbol = []
            for atom in atoms2:
                at.append(atom)
                single_atom = atom.GetSymbol()
                atom_symbol.append(single_atom)
            txt = e.readlines()
            symmetric = 0
            not_symmetric = 0
            linenumber = 0
            up_atom = []
            down_atom = []
            length = len(at)-1
            while linenumber in range(length):
                a = str(txt[linenumber + 4])
                aa = a.split("   ")
                aaa = aa[1]
                a1 = float(aaa)
                b = str(txt[linenumber + 5]) 
                bb = b.split("   ")
                bbb = bb[1]
                b1 = float(bbb)
                if a1 == b1:
                    symmetric = symmetric + 1
                    up_atom.append(str(linenumber))
                    down_atom.append(str(linenumber + 1))
                    linenumber = linenumber + 2
                else:
                    not_symmetric = not_symmetric + 1
                    if not_symmetric >1:
                        break
                    central_atom = str(linenumber)
                    linenumber = linenumber + 1
                    a = str(txt[linenumber + 4])
                    aa = a.split("   ")
                    aaa = aa[1]
                    a1 = float(aaa)
                    b = str(txt[linenumber + 5])
                    bb = b.split("   ")
                    bbb = bb[1]
                    b1 = float(bbb)
                    if a1 == b1:
                        symmetric = symmetric + 1
                        up_atom.append(str(linenumber))
                        down_atom.append(str(linenumber + 1))
                        linenumber = linenumber + 2
                    else:
                        not_symmetric = not_symmetric + 1
            central_atom_total.append(central_atom)
            up_atom_total.append(up_atom)
            down_atom_total.append(down_atom) 
            symmetric_total.append(symmetric)
            not_symmetric_total.append(not_symmetric)
            
            mm = MWAI(mol2) # Draw a diagram showing the atomic index
            
            TotalAtomBondNum = []
            AtomPairs = []
            AtomIdx = []
            shuju = []
            
            if symmetric < 2:
                not_symmetrical_molecule.append(filename)
        # Checking the number of bonds per atom in each of the 30 molecular structures, and screening out the eligible atomic index
                ff = int(ff)
                ff = ff + 1
                jj = mol2.GetAtoms()
        
                for kk in jj:
                  ll = kk.GetDegree()
                  if 1<ll<3:
                    hh = kk.GetIdx()
                    AtomIdx.append(hh)
                
                # Arranging and combining the data
                  for m in itertools.combinations(AtomIdx,2):
                    shuju.append(m)   
                    shuju = list(set(shuju)) # Removing duplicates
                AtomPairs.append(shuju)
                TotalAtomBondNum.append(AtomIdx)
               
                 # Reading the atomic serial number that needs to be highlighted
                for o in shuju:
                  x1 = o[0]
                  x2 = o[1]
                  p = [x1,x2]
                   
                 # Highlighting the atoms
                  ddd = rdMolDraw2D.MolDraw2DCairo(500, 500) # Width in front, length in back
                  rdMolDraw2D.PrepareAndDrawMolecule(ddd, mm, highlightAtoms = list(p))
                  ddd.FinishDrawing()
                  filename_number = os.path.splitext(filename)
                  filename_number_0 = filename_number[0]
                  if int(filename_number_0) < 30:
                      pass
                  else:
                      ee = int(ee)
                      ee = ee + 1
                      ee = str(ee)
                      ddd.WriteDrawingText(best_edges_png + '/' + 'best_mol_' + filename_number_0 + '_' + ee + '.png')
                      
                      # Modifing the CIF file
                      old_cif = open(best_edges_cif + '/' + filename_number_0 + '.cif' , "r")
                      txts3 = old_cif.read()
                      old_cif.close()
                      point_atom1 = atom_symbol[x1] + str(int(x1)+1)
                      new_atom1 = str(str('X')+str(int(x1)+1))
                      point_atom2 = atom_symbol[x2] + str(int(x2)+1)
                      new_atom2 = str(str('X')+str(int(x2)+1))
                      changed = str(txts3).replace(point_atom1,new_atom1).replace(point_atom2,new_atom2)
                      with open(best_edges_cif + '/' + 'best_mol_' + filename_number_0 + '_' + ee + '.cif',"w") as f2:
                          changed.rstrip() 
                          f2.write(changed)
                      f2.close()
    
            else:
                symmetrical_molecule.append(filename)
                ff = int(ff)
                ff = ff + 1
                for m in itertools.combinations(up_atom,2):
                  shuju.append(m)   
                  shuju = list(set(shuju)) # Removing duplicates
                AtomPairs.append(shuju)
                TotalAtomBondNum.append(AtomIdx)
               
                # Reading the atomic serial number that needs to be highlighted
                for o in shuju:
                  x1 = o[0]
                  x1 = int(x1)
                  x2 = o[1]
                  x2 = int(x2)
                  p = [x1,x2]
        
                 # Highlighting the atoms
                  ddd = rdMolDraw2D.MolDraw2DCairo(500, 500) # Width in front, length in back
                  rdMolDraw2D.PrepareAndDrawMolecule(ddd, mm, highlightAtoms = list(p))
                  ddd.FinishDrawing()
                  filename_number = os.path.splitext(filename)
                  filename_number_0 = filename_number[0]
                  if int(filename_number_0) < 30:
                      pass
                  else:
                      ee = int(ee)
                      ee = ee + 1
                      ee = str(ee)
                      ddd.WriteDrawingText(best_edges_png + '/' + 'best_mol_' + filename_number_0 + '_' + ee + '.png')
                      
                # Modifing the CIF file
                      old_cif = open(best_edges_cif + '/' + filename_number_0 + '.cif' , "r")
                      txts3 = old_cif.read()
                      old_cif.close()
                      point_atom1 = atom_symbol[x1] + str(int(x1)+1)
                      new_atom1 = str(str('X')+str(int(x1)+1))
                      point_atom2 = atom_symbol[x2] + str(int(x2)+1)
                      new_atom2 = str(str('X')+str(int(x2)+1))
                      changed = str(txts3).replace(point_atom1,new_atom1).replace(point_atom2,new_atom2)
                      with open(best_edges_cif + '/' + 'best_mol_' + filename_number_0 + '_' + ee + '.cif',"w") as f2:
                           changed.rstrip() 
                           f2.write(changed)
                      f2.close()