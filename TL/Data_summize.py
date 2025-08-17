# -*- coding: utf-8 -*-

import pandas as pd
import json
import os
import numpy as np
import shutil
import sys

# Defining the file path to store the count
counter_file_path = 'counter.json'
output_excel = 'final_data.xlsx'

# Attempts to read the current count value from the file, initializing it to 0 if the file does not exist
try:
    with open(counter_file_path, 'r') as file:
        data = json.load(file)
        counter = data['counter']
except (FileNotFoundError, json.JSONDecodeError):
    counter = 1

# Adding 1 to the count
counter += 1

df1 = pd.read_excel('TL_data_target_task_train.xlsx', engine='openpyxl')
df2 = pd.read_excel('TL_data_target_task_test.xlsx', engine='openpyxl')

combined_df = pd.concat([df1, df2], ignore_index=True)

columns = combined_df.columns.tolist()
if len(columns) >= 2:
    columns[-2] = 'TSN_from_TL'
    columns[-1] = 'TSN_from_RASPA'
    combined_df.columns = columns
    
sheet_name = f'Cycle {counter}'

# Checking whether the file exists
if os.path.exists(output_excel):
    # If the file exists, append a new worksheet
    with pd.ExcelWriter(output_excel, mode='a', engine='openpyxl') as writer:
        combined_df.to_excel(writer, sheet_name=sheet_name, index=False)
else:
    # If the file does not exist, create a new Excel file and write to it
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        combined_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
# Writing the updated count back to the file
with open(counter_file_path, 'w', encoding='utf-8') as file:
    json.dump({'counter': counter}, file)
    

# Extracting the TSN and judging if the loop should be continued
excel_file_path = 'final_data.xlsx'
sheet_name_1 = f'Cycle {counter}'
sheet_name_2 = f'Cycle {counter-1}'

df_last = pd.read_excel(excel_file_path, sheet_name = sheet_name_2, engine='openpyxl')
df_now = pd.read_excel(excel_file_path, sheet_name = sheet_name_1, engine='openpyxl')

name_now = df_now.iloc[:, 0].values    
TSN_now_max = max(df_now['TSN_from_TL'].values)
if sheet_name_2 == 'Cycle 1':
    TSN_last_max = max(df_last['TSN'].values)
else:
    TSN_last_max = max(df_last['TSN_from_TL'].values)

# Judging
if TSN_last_max <= TSN_now_max:
    
    # Extracting the linkers from the best MOFs, and copy these linkers from best_edges_mol to optimal_linker
    new_MOF_index = np.argsort(df_now['TSN_from_TL'].values, axis = 0)[::-1]
    best_MOF_index = new_MOF_index[:50]
    for i in os.listdir('optimal_linker'):
        os.remove('optimal_linker/' + i)
    
    for i in best_MOF_index:
        best_MOF_name = name_now[i]
        best_MOF_linker = best_MOF_name.split('best_mol_')[1].split(' ')[0] + '.mol'
        
        # Returning 'ADTL/' to copy the linkers
        current_directory = os.getcwd()
        parent_directory = os.path.dirname(current_directory)
        os.chdir(parent_directory)
        
        shutil.copy('best_edges_mol/' + best_MOF_linker, 'optimal_linker/' + best_MOF_linker)
        
    # Outputing status code
    sys.exit(0)
    
else:
    sys.exit(1)