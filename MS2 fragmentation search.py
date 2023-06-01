# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:17:30 2023

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
os.chdir("c:/users/mudoe/desktop/hetexp working/h5")
files = glob.glob(os.path.join('*.h5'))
print(len(files))
files = [f for f in files if not '_NEG' in f]
print(len(files))
df_files = pd.DataFrame()
df_files['filename'] = files
df_files['run_order'] = df_files['filename'].apply(lambda x: os.path.basename(x).split('_')[-1].replace('.h5','').replace('Run',''))
df_files['run_order'] = df_files['run_order'].astype(int)
df_files.sort_values('run_order',ascending=True,inplace=True)
df_files.reset_index(inplace=True,drop=True)
df_files

#%%features of interest
os.chdir("c:/users/mudoe/desktop/hetexp working")
qu = pd.read_csv("Aaron query.csv")
qu = qu.where(qu["_dyn_#rt"] < 10)
qu = qu.where(qu["_dyn_#rt"] > .7)
qu = qu.dropna(subset=["rowid"])
qu = qu.reset_index()
qu = qu.drop("index",axis=1)

#%%crunch
###
###this is the block/loop to edit if you want to crunch differently
###
qu2 = copy.deepcopy(qu)
for i in range(len(qu2)):
    for j in range(len(qu2)):
        if i !=j:
            if (qu2["_dyn_#precmz"][i] - (qu2["_dyn_#precmz"][i]/200000)) < qu2["_dyn_#precmz"][j] < (qu2["_dyn_#precmz"][i] + (qu2["_dyn_#precmz"][i]/200000)):
                x = qu2.iat[i,4]-.15
                y = qu2.iat[j,4]
                z = qu2.iat[i,4] + .15
                if x < y < z:
                    qu2.iloc[j,:] = np.nan
qu2 = qu2.dropna(subset = ["rowid"])
qu2 = qu2.reset_index()
qu2 = qu2.drop("index", axis=1)

#%%features of interest
#df_standards = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/337crunched.csv")
df_standards = copy.deepcopy(qu2)
df_standards = df_standards[["rowid","_dyn_#rt","_dyn_#precmz"]] ##this changes depending on input
nastd = copy.deepcopy(df_standards)
nastd["_dyn_#precmz"] = nastd["_dyn_#precmz"]-H+Na
for i in range(len(nastd)):
    nastd.iat[i,0] = "zNa" + str(nastd.iat[i,0])
df_standards = [df_standards,nastd]
df_standards=pd.concat(df_standards)
df_standards = df_standards.reset_index()

ppm_tolerance = 5.0
extra_time = 0

df_standards.rename(columns={'rowid':'label','_dyn_#precmz':'mz','_dyn_#rt':'rt_peak'},inplace=True)
#df_standards.drop(columns=['Unnamed: 0'],inplace=True)
df_standards['rt_min'] = df_standards['rt_peak'] - 0.1
df_standards['rt_max'] = df_standards['rt_peak'] + 0.1
df_standards['ppm_tolerance'] = ppm_tolerance
df_standards['extra_time'] = extra_time

df_standards['group_index'] = ft.group_consecutive(df_standards['mz'].values[:],
                                         stepsize=ppm_tolerance,
                                         do_ppm=True)

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
    
#%% intensities of features
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
os.chdir("c:/users/mudoe/desktop/hetexp working/h5")
results = []
#for i in range(15,45):
for i in range(12,875):
    x = data_list[i]
    d = ft.get_data(x,return_data=True,save_file=False)
    d['ms1_data']['filename'] = x['lcmsrun']
    d['ms1_data']['group'] = os.path.basename(x['lcmsrun']).split('_')[12]
    d['ms1_data']['sample_blank'] = os.path.basename(x['lcmsrun']).split('_')[14]
    d['ms1_data']['run_order'] = x['file_index']
    lead = d["ms2_data"]
    d = d['ms1_data']
    quick ="FPS" in d.at[0,'filename']
    if quick == False:
        sod = d['label'].str.contains('Na')
        sod = sod.fillna(value=False)
        sod2 = d[sod]
        d = d[~sod]
        lead = lead.loc[lead['label'].isin(d['label']),:]
        for i in range(len(lead['label'].unique())):
            lead2 = lead[lead["label"] == lead['label'].unique()[i]]
            if len(lead2) > 0:
                peak = max(lead2["i"])
                ones = lead2[lead2["mz"] > (102.055 - .05)]
                ones = ones[ones["mz"] < (102.055 + .05)]
                ones = ones.sort_values(by = ["i"], ascending=False)
                eights = lead2[lead2["mz"] > (84.0444 - .05)]
                eights = eights[eights["mz"] < (84.0444 + .05)]
                sevs = lead2[lead2["mz"] > (74.06 - .05)]
                sevs = sevs[sevs["mz"] < (74.06 + .05)]
                fives = lead2[lead2["mz"] > (56.05 - .05)]
                fives = fives[fives["mz"] < (56.05 + .05)]
                if len(ones) > 0:
                    if ones.iat[0,3] > (peak * .3):
                        boolcheck = [len(eights) > 0,len(sevs) > 0,len(fives) > 0]
                        if sum(boolcheck) > 1:
                            passer = 0
                        else:
                            lead.loc[lead["label"] == lead['label'].unique()[i],"rt"] = np.nan
                    else:
                        lead.loc[lead["label"] == lead['label'].unique()[i],"rt"] = np.nan
                else:
                    lead.loc[lead["label"] == lead['label'].unique()[i],"rt"] = np.nan
        lead = lead.dropna(axis = 0)
        d = d.loc[d["label"].isin(lead["label"]),:]
        d = pd.concat([d,sod2])
    results.append(d)
results = pd.concat(results)
print("--- %s seconds ---" % (time.time() - start_time))

#%%skips
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

#results = pd.read_csv("c:/users/mudoe/desktop/hetexp working/bleed results/20220606HNa.csv")
#or
results = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/fresh run/resultsms2.csv")
results = results.drop("Unnamed: 0",axis=1)

#%%something
cols = ['label','filename', 'group',
       'sample_blank', 'run_order']
g = results.groupby(cols).agg(max=('i',np.max),counts=('i',len))
g.reset_index(drop=False,inplace=True)
#Na2 = g["label"].str.contains("Na")
#gH = g[~Na2]
#gNa = g[Na2]
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

#%%solve sodiated issue
indus = gpos.index.values
valid = []
gpos["temp"] = indus
for i in range(len(gpos)):
    for j in range(len(gpos)):
        if i != j:
            if "Na" in indus[i]:
                x = indus[i].split("Na")[1]
            else:
                x = indus[i]
            if "Na" in indus[j]:
                y = indus[j].split("Na")[1]
            else:
                y = indus[j]
            if x == y:
                valid.append(indus[i])
                valid.append(indus[j])
valid = pd.DataFrame(valid)
flounder = valid[0].unique()
z = gpos["temp"].isin(flounder)
gpos = gpos.loc[z,:]
gpos = gpos.drop("temp",axis=1)


#%%labels
def getdf(g):
    df = g[["filename","run_order"]]
    dfps = df['filename'].str.contains("FPS")
    dpos = df[~dfps]
    #dall = copy.deepcopy(df)
    df = dpos
    df = df.sort_values("run_order")
    for i in range(1,len(df)):
        if df.iat[i,1] == df.iat[i-1,1]:
            df.iat[i-1,1] = np.nan
    df = df.dropna(subset = ["run_order"])
    for i in range(len(df)):
        if df.iat[i,0].find("FPS") == -1:
            x = df.iat[i,0]
            if x.find("Batch304-") != -1:
                x = x.split("Batch304-",1)[1]
                x = x.split("_",1)[0]
                df.iat[i,0] = x
            if x.find("Std") != -1:
                x = x.split("Std",1)[1]
                x = x.split("_",1)[1]
                x = x.split("_",1)[0]
                df.iat[i,0] = x+" HSL standard"
    df.rename(columns = {'filename':'plasmid'}, inplace = True)
    df["locus tag"]=1
    df["organism"]=1
    user = pd.read_csv("c:/users/mudoe/desktop/hetexp working/user.csv")
    for i in range(len(user)):
        user.iat[i,0] = user.iat[i,0].split("_",1)[1]
    for i in range(len(df)):
        for j in range(len(user)):
            if df.iat[i,0] == user.iat[j,0]:
                df.iat[i,2] = user.iat[j,1]
    locus = pd.read_csv("c:/users/mudoe/desktop/hetexp working/locus.csv")
    for i in range(len(df)):
        for j in range(len(locus)):
            if df.iat[i,2] == locus.iat[j,2]:
                df.iat[i,2] = locus.iat[j,0]
                df.iat[i,3] = locus.iat[j,1]
    return df
df = getdf(g)
#%%snr
def control():
    alp = 1152
    bet = 1159
    gam = 1149
    delt = 1162
    control = gpos[[alp,bet,gam,delt]]
    control = control.max(axis=1)
    return(control)
control = control()
background_ratio = 10
for j in range(len(posrun)):
    for i in range(len(gpos)):
        if gpos.iat[i,j] < (background_ratio*control.iat[i]):
            gpos.iat[i,j] = 0
            
#%%median thresh
# standards = gpos.index.values
# median = []
# for i in range(len(standards)):
#     median.append(np.median(gpos.iloc[i,:]))
# for j in range(len(posrun)):
#     for i in range(len(gpos)):
#         if gpos.iat[i,j] < median[i]:
#             gpos.iat[i,j] = 0

#%%E6 thresh
for i in range(len(gpos)):
    for j in range(len(gpos.columns)):
        if gpos.iat[i,j] < (5*(10**6)):
            gpos.iat[i,j] = 0

#%%Na       
t = int(len(gpos)/2)
for i in range(t):
    for j in range(len(gpos.columns)):
        if (gpos.iat[i,j]*.1) > gpos.iat[i+t,j]:
            gpos.iat[i,j] = 0
            
#%%bleed
cut = []
standards = gpos.index.values
fstd = gfps.index.values
zeta = np.isin(standards,fstd)
lineup = gpos.loc[zeta,:]
test = lineup.index.values
zeda = np.isin(fstd, test)
gfps = gfps.loc[zeda,:]
for j in range(1,len(posrun)):
    for i in range(len(lineup)):
        if lineup[posrun[j]][i] <= gfps[posrun[j]-1][i]*10:
            if lineup[posrun[j]][i] > (5*(10**5)):
                cut.append([standards[i],df.iat[j,2],df.iat[j,3],posrun[j],lineup[posrun[j]][i],gfps[posrun[j]-1][i]])  
            lineup[posrun[j]][i] = 0
cut = pd.DataFrame(cut)
cut.rename(columns={0:'feature',1:'locus tag',2:'organism',3:"run",4:'sample intensity',5:'blank intensity'},inplace=True)

for i in range(len(gpos)):
    for j in range(len(lineup)):
        if standards[i] == test[j]:
            gpos.iloc[i,:] = lineup.iloc[j,:]
            
#%%output
outs = int(len(gpos)/2)
output = pd.DataFrame(index=range(len(gpos.columns)),columns=range((outs*3)+4))
output[0] = df["organism"].tolist()
output[1] = df["locus tag"].tolist()
output[2] = df["plasmid"].tolist()
output[3] = df["run_order"].tolist()
gpos2 = gpos.iloc[range(int(len(gpos)/2)),:]
for i in range(len(gpos2.columns)):
    holder = pd.DataFrame(index=range(int(len(standards)/2)),columns=["std","i","%i"])
    holder["std"] = standards[range(int(len(standards)/2))]
    holder["i"] = gpos2.iloc[:,i].tolist()
    holder = holder.sort_values("i",ascending=False)
    for j in range(len(holder)):
        if holder.iat[0,1] != 0:
            holder.iat[j,2] = holder.iat[j,1]/holder.iat[0,1]
    holder = holder.to_numpy().flatten()
    output.iloc[i,range(4,(outs*3)+4)] = holder
    
#%%skips
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
Na = 22.989218
H = 1.007276
output = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/fresh run/fillerfilterms2_badplasremoved.csv")
output = output.drop("Unnamed: 0",axis=1)

#%%PMI MOB
pm = pd.read_csv("c:/users/mudoe/desktop/hetexp working/moborpmi.csv")
for i in range(len(output)):
    hold = pm.loc[pm["plasmid"] == output.iat[i,5],:]
    if len(hold) > 0:
        output.iat[i,2] = hold.iat[0,1]
        
#%%KO
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
x = pd.read_csv("c:/users/mudoe/desktop/ko numbers.csv")
for i in range(len(x)):
    z = x.at[i,"Plasmid"] == output['5']
    output.loc[z,"1"] = x.at[i,"KO"]
    if pd.isna(x.at[i,"KO"]):
        output.loc[z,"1"] = str(x.at[i,"KO2"])

#%%NDs
for i in range(len(output)):
    if output.iat[i,8] == 0:
        z = output["5"] == output.iat[i,5]
        output.loc[z,"0"] = "ND"
    
#%%match
for i in range(len(output)):
    if output.iat[i,0] != "ND":
        if output.iat[i,4] != "1":
            check = output["5"] == output.iat[i,5]
            check = output.loc[check,]
            if check.iat[0,7] == check.iat[1,7]:
                output.iat[i,0] = "match"
            else:
                output.iat[i,0] = "inconsistent"

monly = output["0"] == "match"
monly = output.loc[monly,]

#%%skips
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
os.chdir("c:/users/mudoe/desktop/hetexp working")
qu2 = pd.read_csv("c:/users/mudoe/desktop/hetexp working/crunch/2crunch.csv")
qu2 = qu2.drop("Unnamed: 0",axis=1)
qu2 = qu2.where(qu2["_dyn_#rt"] < 10)
qu2 = qu2.where(qu2["_dyn_#rt"] > .7)
qu2 = qu2.dropna(subset=["rowid"])
qu2 = qu2.reset_index()
qu2 = qu2.drop("index",axis=1)
ref = pd.read_csv("c:/users/mudoe/desktop/hetexp working/crunch/ref2.csv")
ququest = pd.read_csv("c:/users/mudoe/desktop/hetexp working/crunch/2crunch.csv")

#%%ref skip only
ref = pd.read_csv("c:/users/mudoe/desktop/hetexp working/crunch/ref2.csv")

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
                qu2.iat[j,0] = "faster" + ref.iat[i,2]
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
    
#%%matched only
mas = []
menos = []
for i in range(len(monly)):
    hold = qu2.loc[qu2["_dyn_#new_scan"] == monly.iat[i,7]]
    if len(hold) < 1:
        menos.append(i)
    elif len(hold) > 1:
        mas.append(i)
    else:
        monly.iat[i,7] = hold.iat[0,0]
print(str(len(menos)) + " matched majors without canon or derived HSL formulae")

#%%all
mas = []
menos = []
for i in range(len(output)):
    hold = qu2.loc[qu2["_dyn_#new_scan"] == output.iat[i,7]]
    if len(hold) < 1:
        mas.append(i)
    elif len(hold) > 1:
        menos.append(i)
    else:
        output.iat[i,7] = hold.iat[0,0]
        
#%%secondary
mas = []
menos = []
for i in range(len(output)):
    hold = qu2.loc[qu2["_dyn_#new_scan"] == output.iat[i,10]]
    if len(hold) < 1:
        mas.append(i)
    elif len(hold) > 1:
        menos.append(i)
    else:
        output.iat[i,10] = hold.iat[0,0]
        
#%%PMI
pmi = pd.DataFrame()
#get matched stats
acp = 0
coa = 0
n = 0
hold = output.loc[output["0"] == "match",]
hold = hold.loc[hold["2"] == "P",]
for i in range(len(hold)):
    if hold.iat[i,1] == "K22955" or hold.iat[i,1] == "K13061" or hold.iat[i,1] == "K13061(2)" or hold.iat[i,1] == "K13060" or hold.iat[i,1] == "K22955(2)":
        acp = acp + 1
    elif hold.iat[i,1] == "K18096" or hold.iat[i,1] == "18096(2)":
        coa = coa + 1
    else:
        n = n + 1
total = acp + coa + n
pmi["match"] = [acp, coa, n, total]
#unmatched stats
acp = 0
coa = 0
n = 0
hold = output.loc[output["0"] == "inconsistent",]
hold = hold.loc[hold["2"] == "P",]
for i in range(len(hold)):
    if hold.iat[i,1] == "K22955" or hold.iat[i,1] == "K13061" or hold.iat[i,1] == "K13061(2)" or hold.iat[i,1] == "K13060" or hold.iat[i,1] == "K22955(2)":
        acp = acp + 1
    elif hold.iat[i,1] == "K18096" or hold.iat[i,1] == "18096(2)":
        coa = coa + 1
    else:
        n = n + 1
total = acp + coa + n
pmi["inconsistent"] = [acp, coa, n, total]
#ND stats
acp = 0
coa = 0
n = 0
hold = output.loc[output["0"] == "ND",]
hold = hold.loc[hold["2"] == "P",]
for i in range(len(hold)):
    if hold.iat[i,1] == "K22955" or hold.iat[i,1] == "K13061" or hold.iat[i,1] == "K13061(2)" or hold.iat[i,1] == "K13060" or hold.iat[i,1] == "K22955(2)":
        acp = acp + 1
    elif hold.iat[i,1] == "K18096" or hold.iat[i,1] == "18096(2)":
        coa = coa + 1
    else:
        n = n + 1
total = acp + coa + n
pmi["ND"] = [acp, coa, n, total]
#total stats
pmi["total"] = pmi.sum(axis = 1)
pmi = pmi/2

#%%quick check
# def check(output):
#     monly = output["0"] == "match"
#     monly = output.loc[monly,]
#     print(str(len(monly["7"].unique())) + " unique major products")
#     pmonly = monly["2"] == "P"
#     pmonly = monly.loc[pmonly,]
# check(output)

#%%add seq to PMI only
#import os
#import sys
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
#from sklearn.manifold import TSNE
#import glob
#import copy
#want to save pdf fonts? then do this:
#import matplotlib
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
#Na = 22.989218
#H = 1.007276
#os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")
#from metatlas.io import feature_tools as ft
#pmi = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/fresh run/ms2withprods.csv")
#seq = pd.read_csv("C:/Users/mudoe/Desktop/hetexp working/plasmid_to_seqs.csv")
#for i in range(len(pmi)):
#    hold = seq[seq["plasmid"] == pmi.iat[i,8]]
#    if len(hold) == 1:
#        pmi.iat[i,0] = hold.iat[0,1]
#        pmi.iat[i,1] = hold.iat[0,2]
