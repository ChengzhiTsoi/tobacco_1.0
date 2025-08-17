# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

# Defining the path of sheet
excel_file = 'final_data.xlsx'
xls = pd.ExcelFile(excel_file, engine='openpyxl')

# Initialization
combined_data_sum_max = []
combined_data_sum_aver = []
combined_data_sum_number = []
tsn_proportion_sum = []
sheet_names = xls.sheet_names

# Reading each worksheet
for sheet in sheet_names:

    data = pd.read_excel(excel_file, sheet_name=sheet, engine='openpyxl')
    
    if 'Cycle 1' in sheet:
        MOF_tsn = data.iloc[:, -1]
    else:
        MOF_tsn = data.iloc[:, -2]
    MOF_name = data["Name"]

    # The number of MOFs
    MOF_number = len(MOF_name)
    combined_data_sum_number.append(MOF_number)

    # Average TSN of MOFs
    MOF_tsn_AVER = MOF_tsn.mean()
    combined_data_sum_aver.append(MOF_tsn_AVER)

    # Max TSN of MOFs
    MOF_tsn_MAX = MOF_tsn.max()
    combined_data_sum_max.append(MOF_tsn_MAX)

    # Defining different performance intervals
    bins = [0, 1, 2, 3, 4, float('inf')]
    labels = ['TSN<1', '1<TSN<2', '2<TSN<3', '3<TSN<4', '4<TSN']
    # Dividing MOFs into different intervals
    data[sheet] = pd.cut(MOF_tsn, bins=bins, labels=labels, include_lowest=True)
    # Calculating the number of MOFs with performance in different intervals
    tsn_distribution = data[sheet].value_counts().sort_index()
    # Calculating the proportion
    tsn_proportion = tsn_distribution / tsn_distribution.sum()
    tsn_proportion_sum.append(tsn_proportion)

# Creating the proportion of DataFrame stores in each performance range
tsn_proportion_df = pd.DataFrame(tsn_proportion_sum, index=sheet_names)
tsn_proportion_df.columns = labels

# Outputing the figure
batches = sheet_names
percentages = tsn_proportion_df.T.to_dict('list')
df = pd.DataFrame(percentages, index=labels)
avg_tsn = combined_data_sum_aver
max_tsn = combined_data_sum_max
mof_counts = combined_data_sum_number

# Beautify Settings
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 20,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 13
})

# Building the figure
fig, ax1 = plt.subplots(figsize = (9, 6))
bar_width = 0.4

# Drawing a stacked bar chart
df.T.plot(
    kind = 'bar',
    stacked = True,
    ax = ax1,
    color = ['lightyellow', 'yellow', 'orange', 'pink', 'purple'],
    edgecolor = 'black',
    width = bar_width
)

ax1.set_xlabel('Batch of MOFs')
ax1.set_ylabel('Proportion of MOFs', labelpad = 15)
ax1.set_title('Performance Distribution and Statistics of MOFs', pad = 20)
ax1.set_ylim(0, 1)
ax1.set_xticks(np.arange(len(batches)))
ax1.set_xticklabels(batches, rotation = 0)

# Secondary y-axes
indices = np.arange(len(batches))

# First right axis: Avg TSN
ax2 = ax1.twinx()
line_avg_tsn, = ax2.plot(indices, avg_tsn, color = 'crimson', marker = 'o', markersize = 8, linewidth = 2, label = 'Avg TSN', zorder = 2)
ax2.set_ylabel('Average TSN', labelpad = 8, color = 'crimson', fontsize = 16)
ax2.tick_params(axis = 'y', colors = 'crimson', pad = 5, labelsize = 14)
ax2.spines['right'].set_color('crimson')
ax2.set_ylim(min(avg_tsn)*0.9, max(avg_tsn)*1.1)

# Second right axis: Max TSN
ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 80))
line_max_tsn, = ax3.plot(indices, max_tsn, color = 'royalblue', marker = 's', markersize = 8, linewidth = 2, label = 'Max TSN', zorder = 2)
ax3.set_ylabel('Max TSN', labelpad = 8, color = 'royalblue', fontsize = 16)
ax3.tick_params(axis = 'y', colors = 'royalblue', pad = 5, labelsize = 14)
ax3.spines['right'].set_color('royalblue')
ax3.set_ylim(min(max_tsn)*0.9, max(max_tsn)*1.1)

# Third right axis: MOF Count
ax4 = ax1.twinx()
ax4.spines['right'].set_position(('outward', 160))
line_mof_counts, = ax4.plot(indices, mof_counts, color = 'forestgreen', marker = '^', markersize = 8, linewidth = 2, label = 'MOF Count', zorder = 2)
ax4.set_ylabel('MOF Count', labelpad = 8, color = 'forestgreen', fontsize = 16)
ax4.tick_params(axis = 'y', colors = 'forestgreen', pad = 5, labelsize = 14)
ax4.spines['right'].set_color('forestgreen')
ax4.set_ylim(min(mof_counts)*0.9, max(mof_counts)*1.1)
ax4.yaxis.set_major_locator(MaxNLocator(nbins = 5))

# Setting the X-axis scale
ax1.set_xticks(indices)
ax1.set_xticklabels(batches)

# Adding legend
lns1, labs1 = ax1.get_legend_handles_labels()
lns2, labs2 = ax2.get_legend_handles_labels()
lns3, labs3 = ax3.get_legend_handles_labels()
lns4, labs4 = ax4.get_legend_handles_labels()
ax1.legend(
    lns1 + lns2 + lns3 + lns4,
    labs1 + labs2 + labs3 + labs4,
    loc = 'upper center',
    bbox_to_anchor = (0.5, -0.15),
    ncol = 3
)

# Dynamically adjust annotation positions
def annotate_points(ax, x, y, color, dynamic_offset, format_str):
    for i, txt in enumerate(y):
        ax.annotate(
            format_str.format(txt),
            (x[i], y[i]),
            textcoords="offset points",
            xytext=dynamic_offset(i, y[i]),
            ha='center',
            color=color,
            fontsize=9,
            bbox=dict(boxstyle = "round, pad = 0.3", edgecolor = color, facecolor = "white"),
            zorder=10
        )
    
# Dynamic offset function to prevent label overlap
def dynamic_offset(i, y_value, data_label):
    if data_label == 'avg':
        return (-20, 12) if i % 2 == 0 else (20, 12)
    elif data_label == 'max':
        return (-20, -15) if i % 2 == 0 else (20, -15)
    elif data_label == 'count':
        return (0, 15) if y_value > 30000 else (0, -15)
    return (0, 10)

# Add all annotations to ensure labels are on the top layer
annotate_points(ax2, indices, avg_tsn, 'red', lambda i, y: dynamic_offset(i, y, 'avg'), '{:.2f}')
annotate_points(ax3, indices, max_tsn, 'blue', lambda i, y: dynamic_offset(i, y, 'max'), '{:.2f}')
annotate_points(ax4, indices, mof_counts, 'green', lambda i, y: dynamic_offset(i, y, 'count'), '{}')

plt.tight_layout()
plt.savefig("Output.png", bbox_inches='tight')
