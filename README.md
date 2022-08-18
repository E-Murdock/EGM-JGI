# EGM-JGI
EGM JGI scripts (R and Python)

novels git.py scans through individual h5 files to find features

comb feature table.py checks the given feature table

massQL processing.py takes a massQL output and converts it to a useful list of HSLs with up to 2 unsaturations in the chain, this script's output is the input for novels git.py

figure.R is multiple blocks to produce figures for this project. One creates a fasta that can be plugged into various programs (CLAP, clustal, etc.) to create a protein dendrogram. One is for creating a heatmap that has the same order as said dendrogram's tips so as to be combined into one figure. The final block is for taking mzmine data and creating overlaid XICs to compare peak RT.

data folder contains a massQL output, a processed massQL output from massQL processing.py, and an unprocessed output from novels git.py
