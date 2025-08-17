# -*- coding: utf-8 -*-

import pandas as pd
import sys
from sklearn.metrics import r2_score

data = pd.read_excel('TL_data_target_task_test.xlsx', engine='openpyxl')
dataset = data.dropna() # Null value

y_simu = dataset.iloc[:, -1].values
y_pred = dataset.iloc[:, -2].values
print("y_simu = ", y_simu)
print("y_pred = ", y_pred)

r2 = r2_score(y_simu, y_pred)
print("R2:", r2)

# Judging whether the model is accuracy
threshold = float(0.6)
if float(r2) >= threshold:
    sys.exit(0)
else:
    sys.exit(1)