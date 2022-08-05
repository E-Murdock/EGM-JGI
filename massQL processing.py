# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:41:17 2022

@author: mudoe
"""

#%% imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import copy
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
Na = 22.989218
H = 1.007276
os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
from metatlas.io import feature_tools as ft
os.chdir("c:/users/mudoe/desktop/hetexp working")
qu = pd.read_csv("massQL Query.csv")
qu = qu.where(qu["_dyn_#rt"] < 10)
qu = qu.where(qu["_dyn_#rt"] > .7)
qu = qu.dropna(subset=["rowid"])
qu = qu.reset_index()
qu = qu.drop("index",axis=1)

#%%get formulae
def mf_finder(mz):
    import requests
    chemcalcURL = 'https://www.chemcalc.org/chemcalc/em'
    options = {'mfRange': 'C5-100H6-100N1-1O3-4',
               'numberOfResultsOnly': False,
               'typedResult': False,
               'useUnsaturation': True,
               'minUnsaturation': 3,
               'maxUnsaturation': 6,
               'jcampBaseURL': 'http://www.chemcalc.org/service/jcamp/',
               'monoisotopicMass': mz,
               'jcampLink': False,
               # The 'jcamplink' returns a link to a file containing isotopic
               # distribution of retrieved molecular structure.
               'integerUnsaturation': True,
               # Ions/Radicals can have non-integer unsaturation
               'referenceVersion': '2013',
               'massRange': (mz/200000)
#              'minMass': -0.5,
#              'maxMass': 0.5,
              }
    return requests.get(chemcalcURL, options).json()
ind = []
molf = []
unsat = []
for i in range(len(qu)):
    table = pd.DataFrame(mf_finder(qu.iat[i,2]-H)['results'])
    if len(table) == 1:
        ind.append(i)
        molf.append(table.iat[0,1])
        unsat.append(table.iat[0,2])
qu2 = qu.loc[ind]
qu2.insert(3,"formula", molf)
qu2.insert(4,"unsat", unsat)
qu2=qu2.reset_index()
qu2=qu2.drop("index",axis=1)

#%%crunch
for i in range(len(qu2)):
    for j in range(len(qu2)):
        if i !=j:
            if qu2["formula"][j] == qu2["formula"][i]:
                x = qu2.iat[i,6]-.15
                y = qu2.iat[j,6]
                z = qu2.iat[i,6] + .15
                if x < y < z:
                    qu2.iloc[j,:] = np.nan
qu2 = qu2.dropna(subset = ["formula"])
qu2 = qu2.reset_index()
qu2 = qu2.drop("index", axis=1)

#%%name 'em
ref = pd.read_csv("c:/users/mudoe/desktop/hetexp working/ref.csv")
ind = []
molf = []
unsat = []
for i in range(int(len(ref))):
    table = pd.DataFrame(mf_finder(ref.iat[i,2]-H)['results'])
    if len(table) == 1:
        ind.append(i)
        molf.append(table.iat[0,1])
        unsat.append(table.iat[0,2])
ref["formula"] = molf
ref["unsat"] = unsat
ref.insert(0,"temp",ref.iloc[:,0])

#%%labelling
carbons = []
hydrogens = []
oxygens = []
for i in range(len(ref)):
    carbons.append(ref["formula"][i].split("H",1)[0])
    hydrogens.append(ref["formula"][i].split("H",1)[1])
    hydrogens[i] = hydrogens[i].split("N",1)[0]
    oxygens.append(ref["formula"][i].split("O",1)[1])
ref["carbons"] = carbons
ref["hydrogens"] = hydrogens
ref["oxygens"] = oxygens

carbons = []
hydrogens = []
oxygens = []
for i in range(len(qu2)):
    carbons.append(qu2["formula"][i].split("H",1)[0])
    hydrogens.append(qu2["formula"][i].split("H",1)[1])
    hydrogens[i] = hydrogens[i].split("N",1)[0]
    oxygens.append(qu2["formula"][i].split("O",1)[1])
qu2["carbons"] = carbons
qu2["hydrogens"] = hydrogens
qu2["oxygens"] = oxygens

for i in range(len(ref)): #find double bonds
    for j in range(len(qu2)):
        if ref.at[i,"carbons"] == qu2.at[j,"carbons"]:
            if ref.at[i,"oxygens"] == qu2.at[j,"oxygens"]:
                if ref.at[i,"row retention time"] > qu2.at[j,"_dyn_#rt"]:
                    if int(ref.at[i,"hydrogens"]) == (int(qu2.at[j,"hydrogens"]) +2):
                        qu2.iat[j,0] = ref.iat[i,2] + " double"
                    if int(ref.at[i,"hydrogens"]) == (int(qu2.at[j,"hydrogens"]) +4):
                        qu2.iat[j,0] = ref.iat[i,2] + " 2 double"
                if ref.at[i,"row retention time"] < qu2.at[j,"_dyn_#rt"]:
                    if int(ref.at[i,"hydrogens"]) == (int(qu2.at[j,"hydrogens"]) +2):
                        qu2.iat[j,0] = ref.iat[i,2] + " unsat but slower"
                    if int(ref.at[i,"hydrogens"]) == (int(qu2.at[j,"hydrogens"]) +4):
                        qu2.iat[j,0] = ref.iat[i,2] + " 2 unsat but slower"

for i in range(len(ref)): #faster, slower, or same as standard (crunched)
#essentially, the faster/slower may be indicative of bad ranges in crunch, branching
#or some that this thing is not an HSL
    for j in range(len(qu2)):
        if ref.iat[i,5] == qu2.iat[j,3]:
            cy = ref.iat[i,4]
            if qu2.iat[j,6] < (cy-.1):
                qu2.iat[j,0] = "faster " + ref.iat[i,2]
            if (cy-.1) < qu2.iat[j,6] < (cy+.1):
                qu2.iat[j,0] = ref.iat[i,2]
            if qu2.iat[j,6] < (cy-.1):
                if ref.iat[i,2].find("oxo") != -1:
                    qu2.iat[j,0] = "potential " + ref.iat[i,2] + " isomer"
            if qu2.iat[j,6] > (cy+.1):
                #if ref.iat[i,2].find("oxo") != -1:
                    #if ref.iat[i,6] == 3:
                        qu2.iat[j,0] = "slower " + ref.iat[i,2]
for i in range(len(qu2)): #add the RT
    if isinstance(qu2.iat[i,0],float):
        qu2.iat[i,0] = qu2.iat[i,3]
    qu2.iat[i,0] = qu2.iat[i,0] + " " + str(qu2.iat[i,6])
###            
###qu2 is your output file, to be fed into novel.git as df_standards
###
