# Import necessary libraries
import pandas as pd
import numpy as np
import seaborn as sns
import math
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

# Setting properties and style to be used by matplotlib 
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica']
rcParams['font.size'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['pdf.fonttype'] = 42

##########################################################################################
# Counts file should be .csv file
#		columns: experimental identifier
#		rows: Fab ID (id)
#		Each cell should contain the counts for a specific FAB ID for a particular 
#		experimental identifier
##########################################################################################
file = '20190829_NGS_matched_counts.csv'

##########################################################################################
# The following exp_dict should be replaced such that:
#		key: value
#		experimental identifier: experimental condition
##########################################################################################
exp_dict = {
	'P3F7': 'Donor 1 \n4h', 
	'P3E8': 'Donor 1 \n4h + \nNK', 'P3D9': 'Donor 1 \n4h + \nIL2 NK',
	'P3C10': 'Donor 1 \n24h', 
	'P3B11': 'Donor 1 \n24h + \nNK', 'P3A12': 'Donor 1 \n24h + \nIL2 NK',
	'P4B2': 'Donor 2 \n4h', 
	'P4C3': 'Donor 2 \n4h + \nNK', 'P4D4': 'Donor 2 \n4h + \nIL2 NK',
	'P4E5': 'Donor 2 \n24h', 
	'P4F6': 'Donor 2 \n24h + \nNK', 'P4G7': 'Donor 2 \n24h + \nIL2 NK'}

##########################################################################################
# Read counts file (set index to id)
##########################################################################################
df = pd.read_csv(file, sep = ',').set_index('id')
#to_drop = [x for x in df.columns if x not in exp_dict]
#df.drop(columns = to_drop, inplace = True)

##########################################################################################
# Normalize to total counts in each condition
# Get log2 fold change between experimental condition and control condition
# Store log2 fold change in norm_df and drop control condition
##########################################################################################
control_dict = {}
control_names = {}
for exp in exp_dict:
	if 'NK' not in exp_dict[exp]:					# Get control condition
		ctrl_df = df[exp].copy()
		ctrl_df_sum = ctrl_df.sum()
		ctrl_df = ctrl_df*100.0/ctrl_df_sum			# Norm to total counts in control condition
		control_dict[exp_dict[exp]] = ctrl_df		# Store normalized control values
		control_names[exp_dict[exp]] = exp			# Store name of control condition

norm_df = df.copy()
for exp in exp_dict:
	ctrl_name = exp_dict[exp].split(' + ')[0]		# Get control to be used
	ctrl_df = control_dict[ctrl_name]
	exp_sum = norm_df[exp].sum()					# Get experimental condition
	norm_df[exp] = norm_df[exp]*100.0/exp_sum		# Norm to total counts in experimental condition
	norm_df[exp] = norm_df[exp]/ctrl_df				# Get fold change as compared to control
	norm_df[exp] = np.log2(norm_df[exp])			# Get log2 fold change

norm_df = norm_df.drop(control_names.values(), axis = 1)
norm_df = norm_df.dropna()
norm_df.rename(columns = exp_dict, inplace = True)

# Reorder experimental groups such that earlier timepoints appear first
don1_names = [col for col in norm_df.columns if 'Donor 1' in col]
norm_don1_df = norm_df[don1_names].sort_index(axis = 1, ascending = False)

don2_names = [col for col in norm_df.columns if 'Donor 2' in col]
norm_don2_df = norm_df[don2_names].sort_index(axis = 1, ascending = False)

norm_df = pd.concat([norm_don1_df, norm_don2_df], axis = 1).reset_index()

##########################################################################################
# Plot heatmap
##########################################################################################
f = plt.figure(figsize=(6.6, 8))
sns.set_style('white')
gs = f.add_gridspec(ncols = 2, nrows = 2, width_ratios = [1, 0.05], height_ratios = [0.3, 1])
ax = f.add_subplot(gs[:, 0])
cbar_ax = f.add_subplot(gs[0, 1])
norm_df = norm_df.set_index('id')
heatmap = sns.heatmap(norm_df, ax = ax, vmin = -2.0, vmax = 2.0, cmap = plt.cm.bwr, 
	cbar_kws = {'ticks':[-2, -1, 0, 1, 2], 'label': 'log' + r'$_2$' + '(fold change)'},
	cbar_ax = cbar_ax,
	xticklabels = True, yticklabels = True,)
	
cbar_ax.tick_params(labelsize=8)
cbar_ax.set_ylabel('log' + r'$_2$' + '(fold change)', size = 8)

# Color yticklabels by donor
for tick_label in ax.get_yticklabels():
	tick_text = tick_label.get_text()
	tick_label.set_fontfamily('sans-serif')
	tick_label.set_fontsize(8)
	tick_label.set_va('center')
	tick_label.set_ha('right')
	if 'CD16' in tick_text:
		tick_label.set_color('c')
	elif 'NCR1' in tick_text:
		tick_label.set_color('teal')
	elif 'NCR3' in tick_text:
		tick_label.set_color('seagreen')
	elif 'TNFRSF9' in tick_text:
		tick_label.set_color('grey')
	elif 'CD244' in tick_text:
		tick_label.set_color('dimgrey')
	if 'GFP' in tick_text:
		tick_label.set_color('g')

# Color xticklabels by antigen
for tick_label in ax.get_xticklabels():
	tick_text = tick_label.get_text()
	tick_label.set_fontfamily('sans-serif')
	tick_label.set_fontsize(8)
	tick_label.set_rotation(0)
	tick_label.set_va('top')
	if 'Donor 1' in tick_text:
		tick_label.set_color('grey')
ax.set_ylabel('')
ax.set_xlabel('')
plt.tight_layout()
plt.savefig('heat.pdf')			# Save heatmap

##########################################################################################
# Correlation plot to show how well experimental conditions correlate
##########################################################################################
def corr_map(x, y, size):
    # Get values such that we will end up with a triangle
    counter = 0
    ind_list = []
    num_of_comparisons = len(x)
    num_of_expts = int(math.sqrt(num_of_comparisons))
    for i in range(num_of_expts):
    	for j in range(num_of_expts):
    		if j <= i:
    			ind_list.append(i*num_of_expts + j)
    		else:
    			break
    x = x.loc[ind_list]
    y = y.loc[ind_list]
    size = size.loc[ind_list]
    
    # Mapping from column names to integer coordinates
    x_labels = [v for v in x.unique()]
    y_labels = [v for v in y.unique()][::-1]

    x_to_num = {p[1]:p[0] for p in enumerate(x_labels)} 
    print(x_to_num)
    y_to_num = {p[1]:p[0] for p in enumerate(y_labels)} 
    print(y_to_num)
    
    size_scale = 500

    f2 = plt.figure(figsize=(6, 6))
    sns.set_style('whitegrid')
    gs2 = f2.add_gridspec(ncols = 1, nrows = 15, hspace=0.2, wspace=0.1) # Setup a 1x15 grid
    ax2 = f2.add_subplot(gs2[:, :])
    
    g = ax2.scatter(x=x.map(x_to_num), # Use mapping for x
    		y=y.map(y_to_num), # Use mapping for y
    		s=size.abs() * size_scale, # Vector of square sizes, proportional to size parameter
    		c=size, # Vector of sizes that are to be mapped to cmap
    		marker='s', # Use square as scatterplot marker
    		cmap = plt.cm.bwr)
    g.set_clim([-1, 1])

    # ...
    # Show column labels on the axes
    ax2.set_xticks([x_to_num[v] for v in x_labels])
    ax2.set_xticklabels(x_labels, rotation = 90, fontsize = 8)
    ax2.set_yticks([y_to_num[v] for v in y_labels])
    ax2.set_yticklabels(y_labels, rotation = 0, fontsize = 8)
    
    ax2.grid(False, 'major')
    ax2.grid(True, 'minor')
    ax2.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax2.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
    
    ax2.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax2.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    
    # Add color legend
    cbar = plt.colorbar(g, orientation = 'vertical', ticks = [-1, -0.5, 0, 0.5, 1])
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_xlabel('Pearson r', size = 12)
    cbar.ax.xaxis.set_label_position('top')
    plt.tight_layout()
    plt.savefig('heat_corr.pdf')			# Save correlation plot

corr = norm_df.corr(method = 'pearson')
corr = pd.melt(corr.reset_index(), id_vars='index') # Unpivot the dataframe, so we can get pair of arrays for x and y
corr.columns = ['x', 'y', 'value']
corr_map(
    x=corr['x'],
    y=corr['y'],
    size=corr['value'])