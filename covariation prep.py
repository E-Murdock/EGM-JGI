# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 10:34:34 2022

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
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#Na = 22.989218
#H = 1.007276
#os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
#rom metatlas.io import feature_tools as ft
import Bio as bio
import Bio.Seq as bs
from Bio import pairwise2
het = pd.read_csv("c:/users/mudoe/desktop/hetexp working/Amy&typemodified.csv")
het = het[het["match type"] < 3]
het = het.dropna(subset=["Amy Naming Nomenclature"])
locus = pd.read_csv("c:/users/mudoe/desktop/hetexp working/locus.csv")
q = bs.Seq("MIVQIGRREEFDKKLLGEMHKLRAQVFKERKGWDVSVIDEMEIDGYDALSPYYMLIQEDTPEAQVFGCWRILDTTGPYMLKNTFPELLHGKEAPCSPHIWELSRFAINSGQKGSLGFSDCTLEAMRALARYSLQNDIQTLVTVTTVGVEKMMIRAGLDVSRFGPHLKIGIERAVALRIELNAKTQIALYGGVLVEQRLAVS")
q2 = pd.DataFrame([*"MIVQIGRREEFDKKLLGEMHKLRAQVFKERKGWDVSVIDEMEIDGYDALSPYYMLIQEDTPEAQVFGCWRILDTTGPYMLKNTFPELLHGKEAPCSPHIWELSRFAINSGQKGSLGFSDCTLEAMRALARYSLQNDIQTLVTVTTVGVEKMMIRAGLDVSRFGPHLKIGIERAVALRIELNAKTQIALYGGVLVEQRLAVS"])

#%%philosophical note
"loop through to do pairwise of each against LasI, we don't care how the seqs line up with each other"
"we  only care how they ine up with LasI (or whatever the query is)"
###query is of course LasI in this case, I guess you could repurpose this if you wanted to

#%%loop through
##parameters
k = (1,-.5,-1,-.2)
def full (het,locus,q):
    locs = het["locus tag"].unique()
    guys = []
    for i in range(len(locs)):
        guys.append([locs[i],locus[locus["Locus Tag"] == locs[i]].iat[0,3],het[het["locus tag"] == locs[i]].iat[0,9]])
    lined = []
###lines below are very important for alignment parameters
    for i in range(len(guys)):
        lined.append([bio.pairwise2.align.globalms(q,bs.Seq(guys[i][1]), k[0], k[1], k[2], k[3])[0][1],bio.pairwise2.align.globalms(q,bs.Seq(guys[i][1]), k[0], k[1], k[2], k[3])[0][0]])
    guys.append(["query","MIVQIGRREEFDKKLLGEMHKLRAQVFKERKGWDVSVIDEMEIDGYDALSPYYMLIQEDTPEAQVFGCWRILDTTGPYMLKNTFPELLHGKEAPCSPHIWELSRFAINSGQKGSLGFSDCTLEAMRALARYSLQNDIQTLVTVTTVGVEKMMIRAGLDVSRFGPHLKIGIERAVALRIELNAKTQIALYGGVLVEQRLAVS","3OC12(lasI)"])
    lined.append(["MIVQIGRREEFDKKLLGEMHKLRAQVFKERKGWDVSVIDEMEIDGYDALSPYYMLIQEDTPEAQVFGCWRILDTTGPYMLKNTFPELLHGKEAPCSPHIWELSRFAINSGQKGSLGFSDCTLEAMRALARYSLQNDIQTLVTVTTVGVEKMMIRAGLDVSRFGPHLKIGIERAVALRIELNAKTQIALYGGVLVEQRLAVS","blank"])
    guys = pd.DataFrame(guys)
    lads = pd.DataFrame(lined)
    full = pd.concat([guys, lads], axis=1)
    full.columns = ["locus","unaligned","signal","aligned","query aligned"]
    return full
full = full(het,locus,q)

#%%find the guys of interest
las = [141,144,156,101,139,102,151,143,19,68,45,98,180,158,74,124]
las.sort()
las2 = []
for i in range(len(las)):
    las2.append(las[i]+1)
sleet = []
for i in range(len(full)-1):
    x = [*full.iat[i,3]]
    y = [*full.iat[i,4]]
    n=0
    pertains = []
    lastemp = pd.DataFrame(las2)
    for j in range(len(x)):
        if y[j] != "-":
            n = n+1
        if n in lastemp[0].tolist():
            pertains.append(x[j])
            lastemp[lastemp == n] = "gone"
    sleet.append(pertains)
for i in range(len(sleet)):
    sleet[i] = "".join(sleet[i])
    
#%%see 'em
for i in range(len(sleet)):
    print(sleet[i])

#%%write a fasta?
x = q2.loc[las,:]
x2 = "".join(x[0])
sleet.append(x2) ###get the actual las residues? I guess, that'd probably be worth looking at. what a hassle
full["pertinent"] = sleet
for i in range(len(full)):
    full.iat[i,2] = full.iat[i,2].replace(" ","")

os.chdir("c:/users/mudoe/desktop/covard")
ofile = open("./" + str(k) +".fasta",  "w")
for i in range(len(full)-1):
    ofile.write(">" + str(i)+full.iat[i,2] + "\n" +full.iat[i,5] + "\n")
ofile.close()

#%%get pertinents from other sequences
# jest = "MFVIIQAHEYQKYAAVLDQMFRLRKKVFADTLCWDVPVIGPYERDSYDSLAPAYLVWCNDSRTRLYGGMRLMPTTGPTLLYDVFRETFPDAADLIAPGIWEGTRMCIDEEAIAKDFPEIDAGRAFSMMLLALCECALDHGIHTMISNYEPYLKRVYKRAGAEVEELGRADGYGKYPVCCGAFEVSDRVLRKMRAALGLTLPLYVRHVPARSVVTQFLEMAA"
# j2 = bs.Seq(jest)
# x=bio.pairwise2.align.globalms(q,j2, k[0], k[1], k[2], k[3])[0][1]
# y=bio.pairwise2.align.globalms(q,j2, k[0], k[1], k[2], k[3])[0][0]
# n = 0
# gems = []
# lastemp = pd.DataFrame(las2)
# for j in range(len(x)):
#     if y[j] != "-":
#         n = n+1
#     if n in lastemp[0].tolist():
#         gems.append(x[j])
#         lastemp[lastemp == n] = "gone"
# gems = "".join(gems)

#%%the code below isn't really what we need, but it's kind of interesting and a good reminder of syntax
#temp= []
#for i in range(len(alignment)):
#    temp.append([*alignment.iat[i,2]])
#temp = pd.DataFrame(temp).transpose()

#%%tree stuff
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import Bio.Phylo as Phylo
#zep = AlignIO.read("c:/users/mudoe/desktop/covard/" + str(k) + ".fasta","fasta")
zep = AlignIO.read("c:/users/mudoe/desktop/covard/c6test.fasta","fasta")
calculator = DistanceCalculator('blosum62')
constructor = DistanceTreeConstructor()
dm = calculator.get_distance(zep)
njtree = constructor.nj(dm)
#print(njtree)
#Phylo.draw(njtree)
#Phylo.write(njtree,"c:/users/mudoe/desktop/covard/" + str(k) + ".tre","newick")
Phylo.write(njtree,"c:/users/mudoe/desktop/covard/c6test2.tre","newick")


