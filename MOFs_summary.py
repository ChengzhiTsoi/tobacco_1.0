# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:03:33 2024

@author: Coretshz
"""

import shutil
import os

cycle_number=2
os.makedirs('New_MOF_summary/' + 'Cycle ' + str(cycle_number))

for i in os.listdir('new_designed_mofs'):
    if '.cif' not in i:
        pass
    else:
        shutil.copy('new_designed_mofs/' + i, 'New_MOF_summary/Cycle ' + str(cycle_number) + '/' + i)