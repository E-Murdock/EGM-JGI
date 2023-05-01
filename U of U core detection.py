# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 12:13:16 2022

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
Na = 22.989218
H = 1.007276
os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
from metatlas.io import feature_tools as ft

#%%h5 files
os.chdir("c:/users/mudoe/desktop/Ethan_1209/h5")
files = glob.glob(os.path.join('*.h5'))
files = [f for f in files if not '_NEG' in f]
df_files = pd.DataFrame()
df_files['filename'] = files
df_files['run_order'] = df_files['filename'].apply(lambda x: os.path.basename(x).split('08_')[-1].replace('.h5','').replace('Run',''))
df_files['run_order'] = df_files['run_order'].astype(int)
df_files.sort_values('run_order',ascending=True,inplace=True)
df_files.reset_index(inplace=True,drop=True)

#%%features of interest
os.chdir("c:/users/mudoe/desktop/hetexp working/")
df_standards = pd.read_csv('EGM HSL Core 2.csv')
df_standards = df_standards[["rowid","_dyn_#rt","_dyn_#precmz"]] ##this changes depending on input
nastd = copy.deepcopy(df_standards)
nastd["_dyn_#precmz"] = nastd["_dyn_#precmz"]-H+Na
nastd["rowid"] = "zNa "+nastd["rowid"]
df_standards = [df_standards,nastd]
df_standards=pd.concat(df_standards)
df_standards = df_standards.reset_index()

ppm_tolerance = 5.0
extra_time = 0

df_standards.rename(columns={'rowid':'label','_dyn_#precmz':'mz','_dyn_#rt':'rt_peak'},inplace=True)
#df_standards.drop(columns=['Unnamed: 0'],inplace=True)
df_standards['rt_min'] = df_standards['rt_peak'] - .1
df_standards['rt_max'] = df_standards['rt_peak'] + .1
df_standards['ppm_tolerance'] = ppm_tolerance
df_standards['extra_time'] = extra_time

df_standards['group_index'] = ft.group_consecutive(df_standards['mz'].values[:],
                                         stepsize=ppm_tolerance,
                                         do_ppm=True)
#df_standards.head()
#Nastd = df_standards[df_standards['label'].str.contains("Na")]

#%%data list
data_list = []
for i,row in df_files.iterrows():
    data_setup = {}
    data_setup['lcmsrun'] = row['filename']
    # data_setup['ppm_tolerance'] = ppm_tolerance
    # data_setup['extra_time'] = 0
    data_setup['atlas'] = df_standards
    data_setup['file_index'] = int(os.path.basename(row['filename']).split('_')[-1].replace('.h5','').replace('Run',''))
    data_setup['polarity'] = "positive"
    data_list.append(data_setup)

#%%intensities of features
def unmap_vars(x):
    d = ft.get_data(x,return_data=True,save_file=False)
    d['ms1_data']['filename'] = x['lcmsrun']
    d['ms1_data']['filename'] = os.path.basename(x['lcmsrun'])
    d['ms1_data']['group'] = os.path.basename(x['lcmsrun']).split('_')[12]
    d['ms1_data']['sample_blank'] = os.path.basename(x['lcmsrun']).split('_')[14]
    d['ms1_data']['run_order'] = x['file_index']
    return d['ms1_data']

import time
start_time = time.time()
os.chdir("c:/users/mudoe/desktop/Ethan_1209/h5")
results = []
for i in range(len(data_list)):
    x = data_list[i]
    d = ft.get_data(x,return_data=True,save_file=False)
    d['ms1_data']['filename'] = x['lcmsrun']
    d['ms1_data']['group'] = x['file_index']
    d['ms1_data']['sample_blank'] = x['file_index']
    d['ms1_data']['run_order'] = x['file_index']
    lead = d["ms2_data"]
    d = d['ms1_data']
    results.append(d)
results = pd.concat(results)
print("--- %s seconds ---" % (time.time() - start_time))

#%%skips
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import copy
Na = 22.989218
H = 1.007276
os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
from metatlas.io import feature_tools as ft
results = pd.read_csv("c:/users/mudoe/desktop/Ethan_0117/skips/results.csv")
results = results.drop("Unnamed: 0",axis=1)

#%%something
cols = ['label','filename', 'group',
       'sample_blank', 'run_order']
g = results.groupby(cols).agg(max=('i',np.max),counts=('i',len))
g.reset_index(drop=False,inplace=True)

#%%unmelt
def extractgs(g):
    fps = g['filename'].str.contains("FPS")
    gfps = g[fps]
    gfps = gfps[["label","run_order","max"]]
    fpsrun = gfps["run_order"].unique()
    fpsrun.sort()
    gfps = gfps.pivot_table(index=["label"], columns='run_order')
    gfps.columns = fpsrun
    gfps = gfps.fillna(0)
    gpos = g[~fps]
    gpos = gpos[["label","run_order","max"]]
    posrun = gpos["run_order"].unique()
    posrun.sort()
    gpos = gpos.pivot_table(index=["label"], columns='run_order')
    gpos.columns = posrun
    gpos = gpos.fillna(0)
    runs = g["run_order"].unique()
    runs.sort()
    return (gpos,gfps,posrun)
gpos,gfps,posrun = extractgs(g)
gpospre = copy.deepcopy(gpos)

#%%Na
nagpos = copy.deepcopy(gpos)
t = int(len(gpos)/2)
for i in range(t):
    for j in range(len(gpos.columns)):
        if gpos.iat[i+t,j] < 100:
            nagpos.iat[i,j] = 0
            
#%%output
outs = int(len(gpos)/2)
output = pd.DataFrame(index=range(len(gpos.columns)),columns=range((outs*3)+1))
output[0] = gpos.columns
gpos2 = gpos.iloc[range(int(len(gpos)/2)),:]
standards = gpos.index.values
for i in range(len(gpos2.columns)):
    holder = pd.DataFrame(index=range(int(len(standards)/2)),columns=["std","i","%i"])
    holder["std"] = standards[range(int(len(standards)/2))]
    holder["i"] = gpos2.iloc[:,i].tolist()
    holder = holder.sort_values("i",ascending=False)
    for j in range(len(holder)):
        if holder.iat[0,1] != 0:
            holder.iat[j,2] = holder.iat[j,1]/holder.iat[0,1]
    holder = holder.to_numpy().flatten()
    output.iloc[i,range(1,(outs*3)+1)] = holder

#%%naoutput
outs = int(len(nagpos)/2)
naoutput = pd.DataFrame(index=range(len(nagpos.columns)),columns=range((outs*3)+1))
naoutput[0] = nagpos.columns
nagpos2 = nagpos.iloc[range(int(len(nagpos)/2)),:]
standards = nagpos.index.values
for i in range(len(nagpos2.columns)):
    holder = pd.DataFrame(index=range(int(len(standards)/2)),columns=["std","i","%i"])
    holder["std"] = standards[range(int(len(standards)/2))]
    holder["i"] = nagpos2.iloc[:,i].tolist()
    holder = holder.sort_values("i",ascending=False)
    for j in range(len(holder)):
        if holder.iat[0,1] != 0:
            holder.iat[j,2] = holder.iat[j,1]/holder.iat[0,1]
    holder = holder.to_numpy().flatten()
    naoutput.iloc[i,range(1,(outs*3)+1)] = holder
