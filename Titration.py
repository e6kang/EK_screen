import numpy as np
import math
import re
import os
import colorsys
import pandas as pd
import itertools
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.optimize import curve_fit, leastsq
import seaborn as sns

#########################################################################################
# Dictionaries to be used in functions below
#########################################################################################

IgG_dict = {'27': 'NCR1.11', '72': 'NCR3.12', '78': 'NCR3.18', '3': 'CD16.03',
	'8': 'NCR1.01', '20': 'NCR1.05', '79': 'NCR3.19', '98': 'TNFRSF9.01',
	'108': 'CD16.08'}

Color_dict = {'27': (0, 0, 1), '72': (1, 0, 0), '78': (0.125, 0.25, 1), 
	'3': (0.25, 0.5, 1), '8': (1, 0.2, 0.083), '20': (1, 0.4, 0.167), 
	'79': (0.325, 0.75, 1), '98': (1, 0.6, 0.25), '108': (1, 0.8, 0.33)}

Marker_dict = {'27': ('o', 1), '72': ('o', -1), '78': ('s', 1), '3': ('^', 1),
	'8': ('s', -1), '20': ('^', -1), '79': ('v', 1), '98': ('v', -1),
	'108': ('d', -1)}	

#########################################################################################
# Get file path
#########################################################################################
def getPath(root = "", foldername = ""):
	if root == " " or root == "":
		root = os.getcwd()
	path = root + "/" + foldername
	return path

#########################################################################################
# Create a df based on a .csv file
#########################################################################################
def readCSV (filepath):
	df = pd.read_csv(filepath, lineterminator='\r', error_bad_lines = False)
	Ab_list = ['3', '27', '78', '79', '108', '8', '20', '72', '98']
	Ab_set = set(Ab_list)
	df_cols = df.columns
	
	# Resort columns to match order in Ab_list
	to_sort = set(df_cols).intersection(Ab_list)
	sorted = []
	for i in Ab_list:
		if i in to_sort:
			sorted.append(i)
	
	reordered_cols = []
	count = 0
	for i in df_cols:
		if i not in sorted:
			reordered_cols.append(i)
		else:
			reordered_cols.append(sorted[count])
			count += 1
	
	return df[reordered_cols]

##########################################################################################
# Plot cytotoxicity curves for redirected lysis assays (with IgG)
##########################################################################################
def plot_cytotox(df, color_dict, marker_dict, label_dict):
	rcParams['font.family'] = 'sans-serif'
	rcParams['font.sans-serif'] = ['Helvetica']
	rcParams['mathtext.fontset'] = 'custom'
	rcParams['mathtext.rm'] = 'sans'
	rcParams['pdf.fonttype'] = 42
	
	mean_df = df.groupby('ug/mL').mean().reset_index()
	sem_df = df.groupby('ug/mL').sem().reset_index()

	def PL3_func(x, a, b, c):
		return c/(1.0 +(x/b)**(-a))
	
	# Residuals of actual and fit values
	def residuals(parameters, func, y, x):
		err = y - func(x, *parameters)
		return err
	
	def sum_of_squares(ls):
		s = 0
		for i in ls:
			s += abs(i*i)
		
		return s
	
	sns.set(style = "white", palette = "bright")
	f, ax = plt.subplots(figsize = (3.8, 3))
	
	# Get range of concentrations (get more points to make a nice curve)
	plot_range = np.arange(mean_df['ug/mL'].min(), mean_df['ug/mL'].max(), 0.01)
	
	r_counter = 0
	b_counter = 0
	# For each IgG
	for column in mean_df.columns:
		# Determine if Ab is activating (based on functional screen)
		print(column)
		if column != 'ug/mL':
			if column in color_dict:
				color = color_dict[column]
				marker = marker_dict[column][0]
				if marker_dict[column][1] > 0:
					mfc = color
					mec = color
					mew = 0
				else:
					mfc = 'None'
					mec = color
					mew = 2
			else:
				color = 'black'
				marker = 'o'
				mfc = color
				mec = color
				mew = 0

			# Plot means for each concentration
			if column in label_dict:
				label = label_dict[column]
			else:
				label = column
			ax.errorbar(x = mean_df['ug/mL'], y = mean_df[column], ms = 8,
				xerr = None, yerr = sem_df[column],
				mfc = mfc, mec = mec, mew = mew,
				ecolor = color, fmt = marker, label = label)
			
			
			# Get x and y range
			curr_xrange = mean_df['ug/mL'].tolist()
			curr_yrange = mean_df[column].tolist()
			min0 = 0.02
			max0 = max(curr_yrange)
			
			# Generate fit to means
			ssr = np.inf
			popt = [np.nan, np.nan, np.nan]
			while len(curr_xrange) > 5 and np.isinf(ssr):
				try:
					# Initial guess of parameters
					p0 = [0.1, 0.02, max0]
					# Fit to data with 3PL function
					popt, pcov = leastsq(residuals, p0, args = (PL3_func, curr_yrange, curr_xrange))
					ssr = sum_of_squares(residuals(popt, PL3_func, curr_yrange, curr_xrange))
				except RuntimeError:
					curr_xrange = curr_xrange[:-1]
					curr_yrange = curr_yrange[:-1]
			
			if not np.isinf(ssr):
				PL3_residuals = mean_df[column][:] - PL3_func(mean_df['ug/mL'][:], *popt)
				PL3_ssr = sum_of_squares(PL3_residuals)
				func_to_use = PL3_func
				popt_to_use = popt
			
				# Plot curve fit
				ax.plot(plot_range, func_to_use(plot_range, *popt_to_use), 
					color = color, linewidth = 3)
			
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	
	font_dict = {'fontfamily': 'sans-serif', 'fontsize': 8, 'fontstyle': 'normal'}
	
	ax.set_xlim([0.02, mean_df['ug/mL'].max() + 0.1])
	ax.set_xlabel(r'Concentration ($\mu$g/mL)', fontdict = font_dict)
	
	mean_cytotox = mean_df.drop('ug/mL', axis = 1)
	ax.set_ylim([min(mean_cytotox.min().tolist()) - 5, max(mean_cytotox.max().tolist()) + 5])
	ax.tick_params(labelsize = 8)
	ax.set_ylabel('% cytotoxicity', fontdict = font_dict)
	legend = ax.legend(bbox_to_anchor=(1.01, 1), loc = 2, 
		prop = {'family': 'sans-serif', 
		'size': 8, 'style': 'normal'}, 
		labelspacing = 1.3, handletextpad = 0.25, borderaxespad=0.,
		columnspacing = 0.5, frameon = False)
	plt.tight_layout(pad = 0.5)
	#plt.savefig('mIgG_cytotox.pdf')
	plt.show()

##########################################################################################
# Plot on cell titrations
##########################################################################################
def plot_titration(df, ax, cols, color_dict, marker_dict, label_dict):
	rcParams['font.family'] = 'sans-serif'
	rcParams['font.sans-serif'] = ['Helvetica']
	rcParams['mathtext.fontset'] = 'custom'
	rcParams['mathtext.rm'] = 'sans'
	rcParams['pdf.fonttype'] = 42
	
	def PL4_func(x, a, b, c, d):
		return d + (a-d)/(1.0 + (x/c)**b)
	
	# Residuals of actual and fit values
	def residuals(parameters, func, y, x):
		err = y - func(x, *parameters)
		return err
		
	def sum_of_squares(ls):
		s = 0
		for i in ls:
			s += abs(i*i)
		return s
	
	# Get range of concentrations (get more points to make a nice curve)
	plot_range = np.arange(df['nM'].min(), df['nM'].max(), 0.01)
	
	# For each IgG
	min_mfi = -1
	max_mfi = -1
	for column in cols:
		# Determine if Ab is activating (based on functional screen)
		print(column)
		if column != 'ug/mL':
			if min_mfi < 0:
				min_mfi = df[column].min()
				max_mfi = df[column].max()
			else:
				curr_min = df[column].min()
				curr_max = df[column].max()
				if curr_min < min_mfi:
					min_mfi = curr_min
				if curr_max > max_mfi:
					max_mfi = curr_max
				
			if column in color_dict:
				color = color_dict[column]
				marker = marker_dict[column][0]
				if marker_dict[column][1] > 0:
					mfc = color
					mec = color
					mew = 3
				else:
					mfc = 'None'
					mec = color
					mew = 3
			else:
				color = 'black'
				marker = 'o'
				mfc = color
				mec = color
				mew = 3

			# Plot means for each concentration
			if column in label_dict:
				label = label_dict[column]
			else:
				label = column
			ax.scatter(x = df['nM'], y = df[column], s = 30,
				facecolor = mfc, edgecolors = mec, linewidths = mew,
				marker = marker, label = label)
			
			# Shorten yrange if there is hooking
			curr_xrange = df['nM'].tolist()[:]
			curr_yrange = df[column].tolist()[:]
			min0 = 0
			max0 = max(curr_yrange)
			index_max0 = curr_yrange.index(max0)
			
			curr_xrange = curr_xrange[(index_max0):-1]
			curr_yrange = curr_yrange[(index_max0):-1]
			
			# Generate fit to means
			ssr = np.inf
			popt = [np.nan, np.nan, np.nan, np.nan]
			while len(curr_xrange) > 5 and np.isinf(ssr):
				try:
					# Initial guess of parameters
					p0 = [min0, 1, 1, max0]
					# Fit to data with 4PL function
					popt, pcov = leastsq(residuals, p0, args = (PL4_func, curr_yrange, curr_xrange))
					ssr = sum_of_squares(residuals(popt, PL4_func, curr_yrange, curr_xrange))
				except RuntimeWarning:
					curr_xrange = curr_xrange[:-1]
					curr_yrange = curr_yrange[:-1]
			
			if not np.isinf(ssr):
				PL4_residuals = df[column][1:] - PL4_func(df['nM'][1:], *popt)
				PL4_ssr = sum_of_squares(PL4_residuals)
				func_to_use = PL4_func
				popt_to_use = popt
			
				# Plot curve fit
				ax.plot(plot_range, func_to_use(plot_range, *popt_to_use), 
					color = color, linewidth = 3)
				
	font_dict = {'fontfamily': 'sans-serif', 'fontsize': 8, 'fontstyle': 'normal'}
	
	ax.set_xscale('log')
	ax.set_xlim([0.01, 1500])
	ax.set_xlabel(r'Concentration (nM)', fontdict = font_dict)
	
	mean_mfi = df.drop('nM', axis = 1)
	ax.set_ylim([min_mfi - 1500, max_mfi + 15000])
	ax.tick_params(labelsize = 8)
	ax.set_ylabel('MFI', fontdict = font_dict)
	legend = ax.legend(loc = 2, 
		prop = {'family': 'sans-serif', 
		'size': 8, 'style': 'normal'}, 
		labelspacing = 1, handletextpad = 0.25,
		columnspacing = 0.5)

def plot_all_titrations(df, color_dict, marker_dict, label_dict):
	CD16 = ['3', '108']
	NCR1 = ['27', '8', '20']
	NCR3 = ['78', '79', '72']
	TNFRSF9 = ['98']
	
	sns.set(style = "ticks", palette = "bright")
	f, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (7, 5))
	plot_titration(df, ax[0, 0], CD16, color_dict, marker_dict, label_dict)
	plot_titration(df, ax[0, 1], NCR3, color_dict, marker_dict, label_dict)
	plot_titration(df, ax[1, 0], NCR1, color_dict, marker_dict, label_dict)
	plot_titration(df, ax[1, 1], TNFRSF9, color_dict, marker_dict, label_dict)
	plt.tight_layout(pad = 0.5)
	plt.savefig('NK_MFI.pdf')


##########################################################################################
##########################################################################################
# Run the program
#		Uncomment the first block of code to plot cytotoxicity data
# 		Uncomment the second block of code to plot on cell titration data
##########################################################################################
##########################################################################################
file = '20190502.csv'
cytotox_df = readCSV(getPath(foldername = file))
plot_cytotox(cytotox_df, color_dict = Color_dict, marker_dict = Marker_dict, label_dict = IgG_dict)
"""
file = 'MFI.csv'
mfi_df = readCSV(getPath(foldername = file))
plot_all_titrations(mfi_df, color_dict = Color_dict, marker_dict = Marker_dict, label_dict = IgG_dict)
"""