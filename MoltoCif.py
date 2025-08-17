# -*- coding: utf-8 -*-

content = '''data_GZHULINKER
_audit_creation_date              2023-09-11
_audit_creation_method            'LINKER by TSOI'
_symmetry_space_group_name_H-M    'P1'
_symmetry_Int_Tables_number       1
_symmetry_cell_setting            triclinic
loop_
_symmetry_equiv_pos_as_xyz
  x,y,z
_cell_length_a                    10.0000
_cell_length_b                    10.0000
_cell_length_c                    10.0000
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy\n'''
        
content2 = '''loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type\n'''

def MoltoCif(molfile, ciffile):
    
    from rdkit import Chem
    import datetime
    import math

    
    e = open(molfile)
    e1 = e.readlines()
    number_lines_mol = len(e1)
    mol = Chem.MolFromMolFile(molfile, removeHs=False)
    atom =  mol.GetAtoms()
    atomsnumber = len(atom)
    with open(ciffile, "w+") as oldciffile1:
        current_date = str(datetime.datetime.now().date())
        content1 = content.replace('2023-09-11', current_date)
        oldciffile1.write(content1)
        
        
    i = 0
    while i in range(atomsnumber):
        row = e1[i+4]
        x_coords = row.split()[0]
        y_coords = row.split()[1]
        z_coords = row.split()[2]
        element = row.split()[3]
        row_cif_1 = str(element) + str(i+1)
        row_cif_2 = str(element)
        row_cif_3 = round(float(x_coords)/10, 5)
        row_cif_4 = round(float(y_coords)/10, 5)
        row_cif_5 = round(float(z_coords)/10, 5)
        row_ciffile = row_cif_1 + "       " + row_cif_2  + "       "  + str(row_cif_3) + "    " + str(row_cif_4) + "    " + str(row_cif_5) + '   0.00000  Uiso   1.00\n'
        with open(ciffile, "a+") as oldciffile2:
            oldciffile2.write(row_ciffile)
        i = i + 1
            

    with open(ciffile, "a+") as oldciffile3:
        oldciffile3.write(content2)
        e_oldciffile = open(ciffile)
        lines_cif = e_oldciffile.readlines()
        
    with open(molfile, 'r') as oldmolfile:
        lines_mol = oldmolfile.readlines()
    number_lines_mol = len(lines_mol)

    if 'M  CHG' in lines_mol[number_lines_mol - 2]:
            number_lines_mol = number_lines_mol - 1
    else:
        pass

    i_lines_mol = 0
    while i_lines_mol in range(number_lines_mol - 1 - 4 - atomsnumber):
        lines1_mol = lines_mol[i_lines_mol + 4 + atomsnumber]
        lines12_mol = lines1_mol.split()
        lines121_mol = int(lines12_mol[0])
        lines122_mol = int(lines12_mol[1])
        lines123_mol = lines12_mol[2]
        
        lines1_cif = lines_cif[lines121_mol + 23]
        lines12_cif = lines1_cif.split()
        lines121_cif = lines12_cif[0]
        a_cif = lines121_cif

        lines2_cif = lines_cif[lines122_mol + 23]
        lines22_cif = lines2_cif.split()
        lines221_cif = lines22_cif[0]
        b_cif = lines221_cif
        
        coord11 = float(lines_mol[lines121_mol + 3].split()[0])
        coord21 = float(lines_mol[lines122_mol + 3].split()[0])
        coord12 = float(lines_mol[lines121_mol + 3].split()[1])
        coord22 = float(lines_mol[lines122_mol + 3].split()[1])
        coord13 = float(lines_mol[lines121_mol + 3].split()[2])
        coord23 = float(lines_mol[lines122_mol + 3].split()[2])
        c_cif = math.sqrt((coord11-coord21)*(coord11-coord21)
                          +(coord12-coord22)*(coord12-coord22)
                          +(coord13-coord23)*(coord13-coord23))
           
        if int(lines123_mol) == 1:
            d_cif = "S"
        elif int(lines123_mol) == 2:
            d_cif = "D"
        elif int(lines123_mol) == 3:
            d_cif = "T"
        elif int(lines123_mol) == 4:
            d_cif = "A"
        elif int(lines123_mol) == 5:
            d_cif = "S"
        elif int(lines123_mol) == 6:
            d_cif = "S"
        elif int(lines123_mol) == 7:
            d_cif = "D"
        elif int(lines123_mol) == 8:
            d_cif = "A"

        i_lines_mol = i_lines_mol + 1

        with open(ciffile, "a+") as oldciffile4:
            oldciffile4.write(a_cif + "    " + b_cif + "      " + str(round(c_cif, 3)) + "   .     " + d_cif + '\n')

    # Deleting the last blank line
    with open(ciffile, "r") as newciffile:
        lines_newcif = newciffile.readlines()
        
    number_newcif = len(lines_newcif)

    lines_last = lines_newcif[number_newcif - 1].split('\n')[0]
    e_newcif = open(ciffile)
    e_newcif_content = e_newcif.read()
    e_newcif.close()
    new_e_newcif = str(e_newcif_content).replace(lines_newcif[number_newcif - 1], lines_last)

    with open(ciffile, "w") as newnewciffile:
        newnewciffile.write(new_e_newcif)