# EGM-JGI
EGM JGI scripts (R and Python)

novels git.py scans through individual h5 files to find features

comb feature table.py checks the given feature table

massQL processing.py takes a massQL output and converts it to a list of HSLs with up to 2 unsaturations in the chain, this script's output is the input for novels git.py

figure.R is multiple blocks to produce figures for this project. One creates a fasta that can be plugged into various programs (CLAP, clustal, etc.) to create a protein dendrogram. One is for creating a heatmap that has the same order as said dendrogram's tips so as to be combined into one figure. The final block is for taking mzmine data and creating overlaid XICs to compare peak RT.

waters processing.py searches an LC/MS run that has been processed in MZmine using the methods described here - https://pubs.acs.org/doi/10.1021/acschembio.1c00329 - to search for masses and 102 fragments that correspond to standard HSLs. This can take two inputs, a run with a standard cone voltage for parent masses and a run with an increased cone voltage to find 102 fragments. The blocks are entirely separate, however, so it is possible to only do one. However, providing both will likely give better insights. This does also require a .csv of standard HSL masses and retention times, which will vary by machine.

data folder contains a massQL output, a processed massQL output from massQL processing.py, and an unprocessed output from novels git.py. It also contains the user and locus .csv files that connect JGI runs/plamid numbers to internal construct names, gene locus tags, and organisms of origin.
