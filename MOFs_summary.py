# -*- coding: utf-8 -*-

import shutil
import os

cycle_number=2
os.makedirs('New_MOF_summary/' + 'Cycle ' + str(cycle_number))

for i in os.listdir('new_designed_mofs'):
    if '.cif' not in i:
        pass
    else:
        shutil.copy('new_designed_mofs/' + i, 'New_MOF_summary/Cycle ' + str(cycle_number) + '/' + i)
