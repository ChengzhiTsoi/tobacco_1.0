# -*- coding: utf-8 -*-

content_1 = '''XXX
     RDKit          3D
     
 atom_count bond_count  0  0  chiral_flag  0  0  0  0  0999 V2000/n'''

content_2 = '''M  END

/n'''

def CiftoMol(ciffile, molfile):
    ciffile = "gcmc_selected_mofs/bor_sym_3_mc_0_qzw_sym_4_mc_1_best_mol_64996 549.cif"
    molfile = 'output.mol'
    with open(molfile, "w+") as molfile1:
        molfile1.write(content_1)
        
    import CifFile
    import re
    
    e = open(ciffile, 'r')
    content = e.readlines()
    number = len(content)
    i = 0
    atom_number = 0
    bond_number = 0
    while i in range(number-21):
        row = content[i + 21]
        row1 = row.split()
        if len(row1) > 4:
            element = row1[1]
            
            x = row1[2]
            y = row1[3]
            z = row1[4]
            x1 = x * 10
            y1 = y * 10
            z1 = z * 10
            atom_number = atom_number + 1
            
            row_molfile = '   ' + x1 + '   ' + y1 + '   ' + z1 + ' ' + element + '   0  0  0  0  0  0  0  0  0  0  0  0'
            
            with open(molfile, "a+") as molfile2:
                molfile2.write(row_molfile)
            i = i + 1
            
        else:
            atom1 = row1[0]
            
            atom1_1 = re.findall(r'[0-9]+|[A-Z]+[a-z]+', atom1)
            if len(atom1_1) > 1:
                atom1_2 = atom1_1[1]
            else:
                atom1_2 = atom1_1
                
            atom2 = row1[1]
            atom2_1 = re.findall(r'[0-9]+|[A-Z]+[a-z]+', atom2)
            if len(atom2_1) > 1:
                atom2_2 = atom2_1[1]
            else:
                atom2_2 = atom2_1
                
            bond_type = row1[4]
            if bond_type == 'S':
                bond_type_1 = 1
            elif bond_type == 'D':
                bond_type_1 = 2
            elif bond_type == 'T':
                bond_type_1 = 3
            elif bond_type == 'A':
                bond_type_1 = 4
            
            bond_number = bond_number + 1

            row_molfile = '  ' + atom1_2 + '  ' + atom2_2 + '  ' + bond_type_1 + '  ' + '0'
            
            with open(molfile, "a+") as molfile3:
                molfile3.write(row_molfile)
            i = i + 1


    # Seeking chiral center
    from rdkit import Chem
    import numpy as np
    
    def read_cif_file(cif_file):
        cf = CifFile.ReadCif(cif_file)
        block = cf[cf.keys()[0]]
    
        cell_params = [float(block['_cell_length_a']),
                       float(block['_cell_length_b']),
                       float(block['_cell_length_c']),
                       float(block['_cell_angle_alpha']),
                       float(block['_cell_angle_beta']),
                       float(block['_cell_angle_gamma'])]
    
        symbols = list(block['_atom_site_type_symbol'])
        frac_coords = list(zip(block['_atom_site_fract_x'],
                               block['_atom_site_fract_y'],
                               block['_atom_site_fract_z']))
    
        def frac_to_cart(frac_coords, cell_params):
            a, b, c, alpha, beta, gamma = cell_params
            alpha, beta, gamma = np.radians([alpha, beta, gamma])
            volume = a * b * c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))
            cell_matrix = np.array([[a, b * np.cos(gamma), c * np.cos(beta)],
                                    [0, b * np.sin(gamma), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
                                    [0, 0, volume / (a * b * np.sin(gamma))]])
            cart_coords = np.dot(frac_coords, cell_matrix)
            return cart_coords
    
        cartesian_coords = frac_to_cart(np.array(frac_coords, dtype=float), cell_params)
    
        return symbols, cartesian_coords
    
    def create_rdkit_molecule(symbols, cartesian_coords):
        mol = Chem.RWMol()
        atom_indices = []
        for symbol in symbols:
            atom = Chem.Atom(symbol)
            idx = mol.AddAtom(atom)
            atom_indices.append(idx)
    
        conf = Chem.Conformer(len(atom_indices))
        for i, (x, y, z) in enumerate(cartesian_coords):
            conf.SetAtomPosition(i, (x, y, z))
    
        mol.AddConformer(conf)
        return mol
    
    def find_chiral_centers(mol):
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        return chiral_centers
    
    # reading ciffile
    cif_file = ciffile
    symbols, cartesian_coords = read_cif_file(cif_file)
    
    # Building the molecules
    mol = create_rdkit_molecule(symbols, cartesian_coords)
    
    # Seeking chiral center
    chiral_centers = find_chiral_centers(mol)
    if chiral_centers > 0:
        is_chiral = 1
    else:
        is_chiral = 0
    
    content_mol = ciffile.readlines()
    with open(molfile, "w+") as molfile4:
        content_mol_1 = content_mol.replace('atom_count', atom_number)
        content_mol_2 = content_mol_1.replace('bond_count', bond_number)
        content_mol_3 = content_mol_2.replace('chiral_flag', is_chiral)
        molfile4.write(content_mol_3)
        
    with open(molfile, "w+") as molfile5:
            molfile5.write(content_2)



        