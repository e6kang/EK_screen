# EK_screen
Supplementary scripts to analyze .fasta count files
To process .fasta files and get a .csv file containing count data for further analysis, run match_counts.py. To generate heatmaps and correlation plots with the resulting .csv file, run counts_heatmap.py.

# About match_counts.py
Input: 
  reference file (provided as ref_lib.csv)
  folder directory containing .fasta files
Output:
  .csv file containing all counts that were matched to a sequence in the reference file

# About counts_heatmap.py
Input:
  .csv file containing all counts that were matched to a sequence in the reference file
  (will also need to edit the dictionary in the file to map experimental and control samples)
Output:
  .pdf files of:
      heatmap showing log2 fold changes observed between experimental and control samples
      correlation plots of experimental samples

# About Titration.py and Titration_bispecific
These scripts were used to plot titration curves
Follow comments and headers to use

Scripts written in python 3. 
