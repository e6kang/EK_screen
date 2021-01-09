# Import necessary libraries
import numpy as np
import math
import re
import os
import colorsys
import pandas as pd
import itertools
import warnings
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import leastsq
import seaborn as sns

# Setting properties and style to be used by matplotlib 
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Helvetica'
rcParams['font.size'] = 8
rcParams['pdf.fonttype'] = 42

##########################################################################################	
# Get file path
##########################################################################################
def getPath(root = "", foldername = ""):
	if root == " " or root == "":
		root = os.getcwd()
	path = root + "/" + foldername
	return path

##########################################################################################
# Create a df based on a csv file
##########################################################################################
def readCSV (filepath):
	df = pd.read_csv(filepath, lineterminator='\r', error_bad_lines = False)
	df_cols = df.columns
	cols_to_drop = [x for x in df_cols if ('/' in x) or ('+' not in x and 'nM' not in x)]
	if len(df_cols) - len(cols_to_drop) == 1:
		cols_to_drop = [x for x in df_cols if ('/' in x) and ('nM' not in x)]
	df.drop(columns = cols_to_drop, inplace = True)
	
	curr_cols = df.columns
	new_cols = []
	for c in curr_cols:
		TAA = ''
		NKA = ''
		bi_type = ''
		# Get tumor cell target (either CD20 or Her2)
		if 'CD20' in c:
			TAA = 'CD20'
		elif 'Her2' in c or 'HER2' in c:
			TAA = 'HER2'
		# Get NK cell targeting Ab
		if 'scFv 3' in c:
			NKA = 'CD16.03'
		elif 'scFv 27' in c:
			NKA = 'NCR1.11'
		elif 'scFv 72' in c:
			NKA = 'NCR3.12'
		elif 'scFv 79' in c:
			NKA = 'NCR3.19'
		# Get bispecific type (abbreviation used)
		if 'HC' in c and 'HL' in c:
			bi_type = 'A'
		elif 'HC' in c and 'LH' in c:
			bi_type = 'B'
		elif 'LC' in c and 'HL' in c:
			bi_type = 'C'
		elif 'LC' in c and 'LH' in c:
			bi_type = 'D'
		if len(TAA) > 0:
			if len(NKA) > 0 and len(bi_type) > 0:
				new_cols.append('%sx%s_%s'%(TAA, NKA, bi_type))
			else:
				name = c[1:]
				name = name.split(' + ')[0]
				new_cols.append(name)
		else:
			name = c.split(' + ')[0]
			new_cols.append(name)
			
	#new_cols = [re.sub(r' \+ \w+', '', x) for x in curr_cols]
	df.rename(columns = dict(zip(curr_cols, new_cols)), inplace = True)
	return df
		
##########################################################################################
# Plot on cell titration curves
##########################################################################################
def plot_titration(df, ax, catch_warn = True):
	plot_range = np.arange(0.01, df['nM'].max(), 0.01)
	
	# Lighten rgb value
	def tint(rgb, tint_val):
		rgb = rgb + (tint_val * (1.0 - rgb))
		return rgb
	
	# Darken rgb value
	def shade(rgb, shade_val):
		return rgb * shade_val
	
	# Average of a list of tuples
	def avg_tuple(ls):
		return sum(ls)/len(ls)
	
	# Equation form: y = max + (min - max)/(1+ (x/inflection_pt)**(hill_coeff))
	def PL4_func(x, a, b, c, d):
		return d + (a-d)/(1.0 + (x/c)**b)
			
	# Residuals of actual and fit values
	def residuals(parameters, func, y, x):
		err = y - func(x, *parameters)
		return err
	
	# Sum of squares of a list of numbers
	def sum_of_squares(ls):
		s = 0
		for i in ls:
			s += abs(i*i)
		return s
	# Catch runtime warnings
	if catch_warn:
		warnings.filterwarnings(action = 'error', category = RuntimeWarning)
		
	counter = 0
	cols = df.columns.tolist()[1:]
	
	# For each Ab
	for column in cols:
		# If bispecific not added
		if 'x' not in column:
			if ('CD20' not in column) and ('HER2' not in column):
				color = 'grey'
				marker = 'o'
			else:
				if 'IgG' in column:
					color = 'dimgrey'
					marker = 'o'
				else:
					color = 'black'
					marker = 'o'
		else:
			red, green, blue = 0.0, 0.0, 0.0
			marker = 'X'
			
			# The base_color will depend on NK-targeting arm
			if 'CD16.03' in column:
				green = 0.6
			elif 'NCR3.12' in column:
				red = 0.9
				blue = 0.7
			elif 'NCR3.19' in column:
				red = 0.5
				blue = 0.8
			else:
				blue = 1.0
			
			# Tint the base_color depending on Ab format and get marker
			if '_A' in column:
				tint_val = 0.1
				marker = 's'
			elif '_B' in column:
				tint_val = 0.5
				marker = '^'
			elif '_C' in column:
				tint_val = 0.3
				marker = 'v'
			else:
				tint_val = 0.7
				marker = 'd'
			red = tint(red, tint_val)
			green = tint(green, tint_val)
			blue = tint(blue, tint_val)
			
			# Get the base color
			color = (red, green, blue)
			
		# Initiate mec and mfc and ls
		mec = color
		mfc = color
		ls = '-'
		
		nM_ls = df['nM'].tolist()					# list of concentrations used
		MFI_ls = df[column].tolist()				# list of MFI values in column
		num_Nans = df[column].isna().sum()			# count up number of NaNs in column
		first = float(MFI_ls[0])					# First number in column
		second = float(MFI_ls[1])					# Second number in column
			
		# Skip if not enough data was collected for the column
		if num_Nans > 5:
			continue
				
		# Determine if there is hooking	
		if first < 0.5 * second:
			hook = True
		else:
			hook = False
			
		# Get generate fit to means, pick out concentrations and MFI of interest
		# remove last point (always index -1 since we can't have zero values)
		# may choose to remove first point (index 1 since hook may be seen at high concentrations)
		if hook:
			xrange = nM_ls[1:-1]
			yrange = MFI_ls[1:-1]
		else:
			# If some nans in MFI ls, remove them and their respective concentrations
			if num_Nans > 0:
				xrange, yrange = [], []
				for ind in range(len(MFI_ls)):
					MFI = MFI_ls[ind]
					if np.isnan(MFI):
						continue
					else:
						xrange.append(nM_ls[ind])
						yrange.append(MFI_ls[ind])
			# If no hooking and no nans, use entire ls
			else:
				xrange = nM_ls[:-1]
				yrange = MFI_ls[:-1]
		
		min0 = min(yrange)
		max0 = max(yrange)

		try:
			# Initial guess of parameters
			p0 = [min0, 1, 1, max0]
			# Fit to data with 4PL function
			PL4_popt, PL4_pcov = leastsq(residuals, p0, args = (PL4_func, yrange, xrange))
			PL4_ssr = sum_of_squares(residuals(PL4_popt, PL4_func, yrange, xrange))
		except RuntimeError:
			# If a fit could not be found, set param to nan
			# And sum of squared residuals to inf
			PL4_popt = [np.nan, np.nan, np.nan, np.nan]
			PL4_ssr = np.inf
			
		# If a fit was determined, plot it
		Kd_est = ''
		if not np.isinf(PL4_ssr):
			if PL4_popt[2] > 1000:
				Kd_est = '$K_{D,est}$: %.1f $\mu$M'%(PL4_popt[2]/1000.0)
			else:
				Kd_est = '$K_{D,est}$: %.1f nM'%(PL4_popt[2])
			# Plot curve fit
			ax.plot(plot_range, PL4_func(plot_range, *PL4_popt), 
				color = color, 
				ls = ls, linewidth = 5, dash_capstyle = 'round')
		else:
			Kd_est = '$K_{D,est}$: not found'
		print(column, Kd_est, '\n', 'Bmax: %.1f'%PL4_popt[3])
		# Plot actual data points
		label = column
		ax.scatter(x = df['nM'].tolist(), y = df[column].tolist(), 
			label = label, color = color, marker = marker, edgecolors = mec, s = 100,
			zorder = 5)
	
	font_dict = {'fontfamily': 'sans-serif', 'fontsize': 8, 'fontstyle': 'normal'}
	legend = ax.legend(ncol = 1, loc = 2, prop = {'family': 'sans-serif', 
		'size': 8, 'style': 'normal'}, 
		labelspacing = 0.5, handletextpad = 0.25, 
		columnspacing = 1.8)
	plt.setp(legend.get_texts(), va='center')
	ax.set_xscale('log')
	ax.tick_params(labelsize = 8)
	ax.set_xlabel('Concentration (nM)', fontdict = font_dict)
	ax.set_ylabel('MFI', fontdict = font_dict)
	ax.set_xlim([0.005, 10000])
	ax.set_ylim([0, ax.get_ylim()[1]])

##########################################################################################
# Plot cytotoxicity curves
##########################################################################################
def plot_cytotox(df, ax, catch_warn = False):
	mean_df = df.groupby('nM').mean().reset_index()
	sem_df = df.groupby('nM').sem().reset_index()

	# Lighten rgb value
	def tint(rgb, tint_val):
		rgb = rgb + (tint_val * (1.0 - rgb))
		return rgb
	
	# Darken rgb value
	def shade(rgb, shade_val):
		return rgb * shade_val
	
	# Average of a list of tuples
	def avg_tuple(ls):
		return sum(ls)/len(ls)
		
	# 4 parameter logistic
	def PL4_func(x, a, b, c, d):
		return d + (a-d)/(1.0 + (x/c)**b)
	
	# Residuals of actual and fit values
	def residuals(parameters, func, y, x):
		err = y - func(x, *parameters)
		return err
	
	# Sum of squares of a list of numbers
	def sum_of_squares(ls):
		s = 0
		for i in ls:
			s += abs(i*i)
		return s
				
	# Get range of concentrations (get more points to make a nice curve)
	plot_range = np.arange(mean_df['nM'].min(), mean_df['nM'].max(), 0.01)
	
	# Catch runtime warnings
	if catch_warn:
		warnings.filterwarnings(action = 'error', category = RuntimeWarning)
		
	cols = mean_df.columns.tolist()[1:]
	values_df = mean_df.drop('nM', axis = 1)
	max_val = max(values_df.max())
	min_val = min(values_df.min())
	
	# For each Ab
	for column in cols:
		# If bispecific not added
		if 'x' not in column:
			if ('CD20' not in column) and ('HER2' not in column):
				color = 'grey'
				marker = 'o'
			else:
				if 'IgG' in column:
					color = 'dimgrey'
					marker = 'o'
				else:
					color = 'black'
					marker = 'o'
		else:
			red, green, blue = 0.0, 0.0, 0.0
			marker = 'X'
			
			# The base_color will depend on NK-targeting arm
			if 'CD16.03' in column:
				green = 0.6
			elif 'NCR3.12' in column:
				red = 0.9
				blue = 0.7
			elif 'NCR3.19' in column:
				red = 0.5
				blue = 0.8
			else:
				blue = 1.0
			
			# Tint the base_color depending on Ab format and get marker
			if '_A' in column:
				tint_val = 0.1
				marker = 's'
			elif '_B' in column:
				tint_val = 0.5
				marker = '^'
			elif '_C' in column:
				tint_val = 0.3
				marker = 'v'
			else:
				tint_val = 0.7
				marker = 'd'
			red = tint(red, tint_val)
			green = tint(green, tint_val)
			blue = tint(blue, tint_val)
			
			# Get the base color
			color = (red, green, blue)
			
		# Initiate mec and mfc and ls
		mec = color
		mfc = color
		ls = '-'		

		# Generate fit to means (only if bispec or IgG added)
		if 'x' in column or 'IgG' in column:
			xrange = mean_df['nM'].tolist()
			yrange = mean_df[column].tolist()
			min0 = 0
			max0 = max(yrange)
			
			# Shorten yrange if there is hooking
			index_max0 = yrange.index(max0)
			new_xrange = xrange[:(index_max0 + 1)]
			new_yrange = yrange[:(index_max0 + 1)]
			
			curr_xrange = xrange[:]
			curr_yrange = yrange[:]
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
			
			EC50_est = ''
			print(column, popt)
			if not np.isinf(ssr):
				if popt[2] > 1000:
					EC50_est = '$EC_{50,est}$: %.2f $\mu$M'%(popt[2]/1000.0)
				else:
					EC50_est = '$EC_{50,est}$: %.2f nM'%(popt[2])
				func_to_use = PL4_func
				if len(popt) > 4:
					func_to_use = PL5_func
				# Plot curve fit
				ax.plot(plot_range, func_to_use(plot_range, *popt), 
					color = color, ls = ls, linewidth = 3, dash_capstyle = 'round')
			else:
				EC50_est = '$EC_{50,est}$: not found'
		
			print(column, EC50_est)
			print('Emax: %.2f'%popt[3])
			# Plot means and error bars for each concentration
			label = column
			
			ax.errorbar(x = mean_df['nM'], y = mean_df[column], ms = 8,
				xerr = None, yerr = sem_df[column],
				color = color, ecolor = mec, fmt = marker, mfc = mfc, mec = mec,
				label = label)
				
		else:
			# Plot means and error bars for each concentration
			ax.errorbar(x = mean_df['nM'], y = mean_df[column], ms = 8,
				xerr = None, yerr = sem_df[column],
				color = color, ecolor = color, fmt = marker, mfc = mfc, mec = mec,
				label = column)
	
	font_dict = {'fontfamily': 'sans-serif', 'fontsize': 8, 'fontstyle': 'normal'}
	
	# Show legend
	ncols = 1
	if len(mean_df.columns) > 5:
		ncols = 2
	legend = ax.legend(ncol = ncols, loc = 2, prop = {'family': 'sans-serif', 
		'size': 8, 'style': 'normal'}, 
		labelspacing = 0.5, handletextpad = 0.25, 
		columnspacing = 0.5)
	plt.setp(legend.get_texts(), va='center')
	
	ax.set_xlim([mean_df['nM'].min()/2, mean_df['nM'].max()*2])
	
	multiplier = 1.4
	if len(mean_df.columns) > 7:
		multiplier = 1.6
	ax.set_ylim(min_val - 2.5, max_val * multiplier)
	ax.set_xscale('log')
	ax.tick_params(labelsize = 8)
	ax.set_xlabel('Concentration (nM)', fontdict = font_dict)
	ax.set_ylabel('% cytotoxicity', fontdict = font_dict)


##########################################################################################
##########################################################################################
# Run the program
#		Uncomment the first block of code to plot cytotoxicity data
# 		Uncomment the second block of code to plot on cell titration data
##########################################################################################
##########################################################################################
"""				
##########################################################################################
# Files containing cytotoxicity data
##########################################################################################
file_IgG = 'Trispecifics_PBMC.csv'

IgG_df = readCSV(getPath(foldername = file_IgG))

sns.set(style = "ticks", palette = "bright")
f, ax = plt.subplots(figsize = (3.4, 3))
plot_cytotox(IgG_df, ax, catch_warn = False)
plt.savefig('bispec_cytotox.pdf')
"""

"""
##########################################################################################
# Files containing on cell titration data
##########################################################################################
file_3 = 'CD20x3_FlpIn_MFI.csv'

bispec3_df = readCSV(getPath(foldername = file_3))

sns.set(style = "ticks", palette = "bright")
f, ax = plt.subplots(figsize = (3.4, 3))
plot_titration(bispec3_df, ax, catch_warn = False)
plt.savefig('bispec_mfi.pdf')
"""