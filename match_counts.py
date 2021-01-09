# Import necessary libraries
import os
import re
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

##########################################################################################
# count_file_path: Directory where .fasta files are stored
# 			.fastq.gz files were processed with usegalaxy.org as follows
#			-- Removed sequencing artifacts
#			-- Clipped the following adapter sequence
#			------ TACTGGGGTCAAGGAACCCTGGTCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
#			-- FASTQ Masker by quality score (when quality score < 30) with N
#			-- Collapsed sequences (get counts for each unique sequence)
#			------ These collapsed sequences are exported as .fasta files
# ref_lib_path: The .csv file that maps CDR H3 sequence to Fab ID
##########################################################################################
count_file_path = '20190829_NGS/Galaxy_collapsed'
ref_lib_path = 'ref/ref_lib.csv'

##########################################################################################
# Codon to AA dictionary (to be used in functions below)
##########################################################################################
codon_to_AA_dict = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

##########################################################################################
# Parse .fasta files
#		The 'Name' of the .fasta file should follow the convention of the below example:
#		-- 'xxxxx[xxxx_Name].fasta'
# 		For each .fasta file, create a pandas dataframe to store sequence and count data
#		Store these dataframes in a list
# Returns a list of dataframes
##########################################################################################
def parse_count_file(folderPath):
	# Empty list to contain each df collected from each count file
	DFList = []
	
	# For each count file in the folder
	for inputFileName in os.listdir(folderPath):
		# Get the filepath for the count file
		filepath = os.path.join(folderPath, inputFileName)
		# Get exp name
		name = inputFileName.strip('].fasta').split('_')[-1]
		
		count = 0
		count_dict = {}
		# Get counts for each DNA sequence
		for i, line in enumerate(open(filepath)):
			if i % 2 == 0:
				count = int(line.strip().split('-')[-1])
			else:
				seq = line.strip()
				# Only keep sequences with count > 10 and where quality of each bp is good
				if count > 10 and seq.count('N') < 3:
					count_dict[seq] = count
		# Convert dictionary of seq counts into a dataframe, and store in list
		df = pd.DataFrame.from_dict(count_dict, orient = 'index').reset_index()
		df.columns = ['seq', 'counts']
		df['FileName'] = name
		df = df.sort_values(by = ['counts'], ascending = False)
		DFList.append(df)
	return DFList

##########################################################################################
# Parse reference .csv file
# Returns dictionary
#		key: value
#		Sequence: Fab ID
##########################################################################################
def parse_ref_file (filePath):
	# Read ref file as df
	df = pd.read_csv(os.path.join(os.getcwd(), filePath), error_bad_lines = False)
	df['id'] = df['HC3']+'GAC'
	
	unique_ids = df['id'].unique()
	ref_dict = {}
	for id in unique_ids:
		 ref_dict[id] = ', '.join(df['Protein'].loc[df['id'] == id].tolist())
	
	return ref_dict

##########################################################################################
# Match count sequences to reference dictionary
# Store all count data from all .fasta files in a single dataframe and save as .csv file
##########################################################################################
def match_to_ref(count_df_list, ref_dict):
	# Make dictionary to match CDR H3 AA seq to Ab id
	id_to_AAseq_dict = {}
	id_to_seq_dict = {}
	for k, v in ref_dict.items():
		k_list = re.findall('.{1,3}', k[:-3])
		AA_seq = ''.join(map(lambda x: codon_to_AA_dict[x], k_list))
		id_to_AAseq_dict[v] = AA_seq
		id_to_seq_dict[v] = k[:-3]
	
	# Make a match function (allow N in sequence)
	def match_func(seq):   
		if seq in ref_dict:
			return ref_dict[seq]
		
		elif 'N' in seq:
			if seq.count('N') > 1:
				seq = seq[:-3] + 'GAC'
			seq_A = seq.replace('N', 'A')
			if seq_A in ref_dict:
				return ref_dict[seq_A]
			seq_C = seq.replace('N', 'C')
			if seq_C in ref_dict:
				return ref_dict[seq_C]
			seq_T = seq.replace('N', 'T')
			if seq_T in ref_dict:
				return ref_dict[seq_T]
			seq_G = seq.replace('N', 'G')
			if seq_G in ref_dict:
				return ref_dict[seq_G]
	
	df_dict = {}
	seq_df_list = []
	Fab_df_list = []
	# match sequences to known reference ids
	for df in count_df_list:
		name = df['FileName'].iloc[0]
		df['id'] = df['seq'].apply(lambda x: match_func(x))
		total_counts = df['counts'].sum()
		#df = df[pd.isnull(df['id'])]
		df = df.sort_values(by = ['id'], ascending = False)
		seq_df_list.append(df)
		df = df.dropna(axis = 0)
		matched_total_counts = df['counts'].sum()
		print name, total_counts, float(matched_total_counts)/total_counts
		
		# Sum up all counts on sequences that reference the same ids
		grouped_seq_df = df.groupby(['id'])['counts'].sum()
		grouped_seq_df = grouped_seq_df.reset_index()
		print len(grouped_seq_df)
		# Normalize counts of each HC3 to total number of reads for that barcode

		grouped_seq_df['FileName'] = name
		Fab_df_list.append(grouped_seq_df)
		#print grouped_seq_df
	normed_df = pd.concat(Fab_df_list)
	pivot = normed_df.pivot_table(values = 'counts', index = 'id', 
		columns = 'FileName')
	pivot.to_csv('an_eg.csv')
	return df_dict

##########################################################################################
##########################################################################################	
# Run the program
##########################################################################################
##########################################################################################
count_df_list = parse_count_file(count_file_path)
ref_dict = parse_ref_file(ref_lib_path)
matched_counts_dict = match_to_ref(count_df_list, ref_dict)
