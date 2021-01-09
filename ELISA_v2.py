import numpy as np
import networkx as nx
import math
import re
import os
import pandas as pd
import itertools
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rcParams

def plot_ELISA(df, clone_df):
	rcParams['font.family'] = 'sans-serif'
	rcParams['font.sans-serif'] = ['Helvetica']
	rcParams['font.size'] = 8
	rcParams['xtick.labelsize'] = 8
	rcParams['ytick.labelsize'] = 8
	rcParams['legend.fontsize'] = 8
	rcParams['pdf.fonttype'] = 42
	
	sns.set_style('white')
	# Get names of all the sheets in the excel file
	sheet_name_list = df.book.sheet_names()
	
	def get_well(ls):
		c = math.ceil(ls[1]/2.0)
		name = ls[0] + str(int(c))
		return name
	
	num_sheets = len(sheet_name_list)
	f = plt.figure(figsize=(7, 8.8))
	nrows = math.ceil(num_sheets/2)
	gs = f.add_gridspec(ncols = 2, nrows = nrows)
	
	# For each name get the direct and competition signals
	for i in range(len(sheet_name_list)):
		# Get name of the sheet (should refer to antigen name)
		sheet_name = sheet_name_list[i]
		if sheet_name == 'CD16':
			ratio_cutoff = 0.75
			direct_cutoff = 0.25
		elif sheet_name == 'CD244':
			ratio_cutoff = 0.7
			direct_cutoff = 0.5
		else:
			ratio_cutoff = 0.5
			direct_cutoff = 0.45
		
		# Get values from current sheet
		sheet_df = df.parse(i)
		
		# Remove non-plate data
		sheet_df.dropna(axis='index', how='any', inplace = True)
		sheet_df = sheet_df.reset_index()
		sheet_df = sheet_df.iloc[1:, 2:]			# Only keep relevant columns
		sheet_df.columns = [i+1 for i in range(24)]	# Rename columns
		
		# Get direct and competition
		sheet_df_odd = sheet_df.iloc[::2]			# Separate out odd rows
		sheet_df_odd['row'] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
		dir_comp_df = pd.melt(sheet_df_odd, id_vars = 'row')
		# Get corresponding phage well
		dir_comp_df['well'] = dir_comp_df[['row', 'variable']].apply(get_well, axis = 1)
		dir_comp_df['condition'] = np.where(dir_comp_df['variable']%2 == 1, 'direct', 'competition')
		
		# Get Fc and BSA controls
		sheet_df_even = sheet_df.iloc[1::2]			# Separate out even rows
		sheet_df_even['row'] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
		Fc_BSA_df = pd.melt(sheet_df_even, id_vars = 'row')
		# Get corresponding phage well
		Fc_BSA_df['well'] = Fc_BSA_df[['row', 'variable']].apply(get_well, axis = 1)
		Fc_BSA_df['condition'] = np.where(Fc_BSA_df['variable']%2 == 1, 'Fc', 'BSA')

		# Combine direct, competition, Fc, and BSA in one sheet
		combined_df = pd.concat([dir_comp_df, Fc_BSA_df], axis = 0).drop(['row', 'variable'], axis = 1)
		pivot = combined_df.pivot(index='well', columns='condition', values='value')
		
		# Get passing wells
		pivot['comp_ratio'] = pivot['competition']/pivot['direct']		# Calculate competition/direct
		pivot['pass_comp'] = np.where(pivot['comp_ratio'] < ratio_cutoff, True, False)
		pivot['pass_dir'] = np.where(pivot['direct'] > direct_cutoff, True, False)
		pivot['pass_Fc'] = np.where(pivot['direct'] > pivot['Fc'], True, False)
		pivot['pass'] = pivot['pass_comp'] & pivot['pass_dir'] & pivot['pass_Fc']
		pivot['color'] = np.where(pivot['pass'], '#808080', '#CCCCCC')
		pivot = pivot.sort_values(by = 'color', ascending = False)
		pivot.loc[:, 'Fab ID'] = 'unmade'
		
		# Get color of passing wells
		clones = clone_df[clone_df['Antigen'] == sheet_name].copy()
		clones = clones[clones['pBL347'] == 'yes']
		num_unique = len(clones)
		# Colors chosen from https://xkcd.com/color/rgb/
		palette = sns.xkcd_palette(['red', 'peach', 'gold', 'green', 'blue', 'purple',\
		'pink', 'orange', 'pale yellow', 'pale green', 'pale blue', 'lilac',\
		'magenta', 'salmon', 'yellow', 'light green', 'sky blue', 'light purple',\
		'hot pink', 'rose', 'tangerine', 'aqua', 'royal blue', 'indigo',\
		'light pink', 'bluish purple', 'greyish green', 'barney purple', 'tan', 'brown'])

		counter = 0
		for row in clones.iterrows():
			Fab_ID = row[1]['Fab ID']
			ELISA_wells = row[1]['ELISA wells'].split(',')
			palette_color = palette[counter]
			red = int(palette_color[0]*255)
			green = int(palette_color[1]*255)
			blue = int(palette_color[2]*255)
			color = '#%02x%02x%02x' % (red, green, blue)
			for well in ELISA_wells:
				#print(well)
				pivot.loc[well, 'color'] = color
				pivot.loc[well, 'Fab ID'] = Fab_ID
			counter += 1
		unique_passed = pivot[pivot['Fab ID'] != 'unmade'].drop_duplicates('Fab ID')
		unique_passed = unique_passed.sort_values(by = ['Fab ID'], ascending = True)
		
		#fig, ax = plt.subplots()
		ax = f.add_subplot(gs[i//2, i%2])
		# Plot ELISA signal (color in ones that passed)
		g = ax.scatter(x = pivot['comp_ratio'], y = pivot['direct'], 
			marker = 'o',
			c = pivot['color'], 
			ec = ['black' for x in range(len(pivot))], 
			lw = 1,
			s = 50, zorder = 2)
		
		# Create legend
		label_list = []
		patch_list = []
		for index, row in unique_passed.iterrows():
			label_list.append(row['Fab ID'])
			patch_list.append(mpatches.Circle((0,0), radius = 10, fc = row['color']))
		if len(label_list) > 15:
			ax.legend(patch_list, label_list, loc = 1, ncol = 2, 
				prop = {'family': 'sans-serif', 
				'size': 8, 'style': 'normal'}, 
				labelspacing = 0.25, handletextpad = 0.25, 
				columnspacing = 0.25,
				bbox_to_anchor=(1, 1))
		else:
			ax.legend(patch_list, label_list, loc = 1, 
				prop = {'family': 'sans-serif', 
				'size': 8, 'style': 'normal'},
				labelspacing = 0.25, handletextpad = 0.25,
				bbox_to_anchor=(1, 1))
		
		ax.axhline(y = direct_cutoff, color='k', linestyle='--', zorder = 1)
		ax.axvline(x = ratio_cutoff, color='k', linestyle='--', zorder = 1)
		xmin, xmax = ax.get_xlim()
		ymin, ymax = ax.get_ylim()
		print(xmin, xmax, ymin, ymax)
		
		ax.axvspan(xmin = xmin, xmax = ratio_cutoff,
			ymin = (direct_cutoff-ymin)/(ymax - ymin), ymax = 1, 
			facecolor='r', alpha=0.2, zorder = 0)
		
		font_dict = {'fontfamily': 'sans-serif', 'fontsize': 8, 'fontstyle': 'normal'}
		
		ax.set_xlim(xmin, xmax + 0.65)
		ax.set_ylim(ymin, ymax)
		ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		
		if i//2 == (nrows-1):
			ax.set_xlabel('competition/direct binding', fontdict = font_dict)
		if i%2 == 0:
			ax.set_ylabel('direct binding', fontdict = font_dict)
		ax.set_title('Antigen: ' + sheet_name, fontdict = font_dict)
		
		"""
		passed = pivot[pivot['pass'] == True]
		ax = passed[['direct','competition', 'Fc', 'BSA']].plot(kind='bar', legend=True, fontsize=12)
		ax.set_ylabel('OD450')
		ax.legend(title = '', loc = 'upper left')
		"""
	plt.tight_layout(pad = 0.25, h_pad = 0.5, w_pad = 1.5)
	plt.savefig('ELISA.pdf')

#file = '2018-03-19 CD16, NCR1, NCR3, SF9, SF4.xlsx'
file = 'ELISAs.xlsx'
elisa_data = pd.ExcelFile(file)

clone_file = 'Unique_clones.csv'
clone_data = pd.read_csv(clone_file, error_bad_lines = False)

plot_ELISA(elisa_data, clone_data)