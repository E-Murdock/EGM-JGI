# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:54:16 2022

@author: mudoe
"""

#%%imports
import os
#import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
#from sklearn.manifold import TSNE
import glob
import copy
#want to save pdf fonts? then do this:
import matplotlib
from Bio import SeqIO

#%%data
name = "plasmid.fasta comma.csv"
editing = pd.read_csv("c:/users/mudoe/desktop/hetexp working/puri schaefer fastas/temp/" + name,sep="\t",header = None)
editing[2] = editing[0]
locus = pd.read_csv("c:/users/mudoe/desktop/hetexp working/locus.csv")
match = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/matched_KO2.csv")
maj = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/full_massQL_output.csv")
for i in range(len(editing)):
    editing.iat[i,0] = editing.iat[i,0].replace(">","")
    editing.iat[i,0] = editing.iat[i,0].replace(" ","")
for i in range(len(editing)):
    test = maj.loc[maj["2"] == editing.iat[i,0],]
    if len(test) > 0:
        editing.iat[i,2] = editing.iat[i,2] + " " + test.iat[0,1]
for i in range(len(editing)):
    test = match.loc[match["2"] == editing.iat[i,0],]
    if test.iat[0,1] == "K22955" or test.iat[0,1] == "K13061" or test.iat[0,1] == "K13061(2)" or test.iat[0,1] == "K13060" or test.iat[0,1] == "K22955(2)":
        editing.iat[i,2] = editing.iat[i,2] + " plasmid ACP " + test.iat[0,1]
    elif test.iat[0,1] == "K18096" or test.iat[0,1] == "18096(2)":
        editing.iat[i,2] = editing.iat[i,2] + " plasmid CoA " + test.iat[0,1]
    else:
        editing.iat[i,2] = editing.iat[i,2] + " plasmid undetermined"
for i in range(len(editing)):
    test = match.loc[match["2"] == editing.iat[i,0],]
    if test.iat[0,2] == 1:
        editing.iat[i,2] = editing.iat[i,2] + " matched"
    if test.iat[0,2] == 2:
        editing.iat[i,2] = editing.iat[i,2] + " NDboth"
    if test.iat[0,2] == 3:
        editing.iat[i,2] = editing.iat[i,2] + " inconsistent"
    if test.iat[0,2] == 4:
        editing.iat[i,2] = editing.iat[i,2] + " NDone"
for i in range(len(editing)):
    test = match.loc[match["2"] == editing.iat[i,0],]
    editing.iat[i,2] = editing.iat[i,2] + " " + test.iat[0,4]
    if test.iat[0,3] == "MO enrichment":
        editing.iat[i,2] = editing.iat[i,2] + " enrichment metagenome"
    if test.iat[0,3] == "PMI_endosphere":
        editing.iat[i,2] = editing.iat[i,2] + " endo metagenome"
    if test.iat[0,3] == "PMI_rhizosphere":
        editing.iat[i,2] = editing.iat[i,2] + " rhizo metagenome"
    if test.iat[0,3] == "PMI_soil":
        editing.iat[i,2] = editing.iat[i,2] + " dirt metagenome"     
    if test.iat[0,3] != "MO enrichment" and test.iat[0,3] != "PMI_endosphere" and test.iat[0,3] != "PMI_rhizosphere" and test.iat[0,3] != "PMI_soil":
        editing.iat[i,2] = editing.iat[i,2] + " " +  test.iat[0,3]
        
#editing[0] = editing[0] + " PMI"
editing.to_csv("c:/users/mudoe/desktop/hetexp working/puri schaefer fastas/finished/plasmid organism locus.csv")
