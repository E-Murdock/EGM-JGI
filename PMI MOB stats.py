# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 12:20:42 2022

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
###note that this script determines the statistics when plasmids that did not grow properly are removed. This was determined by OD of E. coli carrying said plasmids.
output = pd.read_csv("c:/users/mudoe/desktop/hetexp working/fullql/fresh run/bad plas removed.csv")
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
                
#%%interrogate majors a bit
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
    
#%%back to it
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

#%%run that on secondary of all
mas = []
menos = []
for i in range(len(output)):
    hold = qu2.loc[qu2["_dyn_#new_scan"] == output.iat[i,10]]
    if len(hold) < 1:
        menos.append(i)
    elif len(hold) > 1:
        mas.append(i)
    else:
        output.iat[i,10] = hold.iat[0,0]
print(str(len(menos)) + " secondary features")

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

#%%mob
mob = pd.DataFrame()
#get matched stats
acp = 0
coa = 0
n = 0
hold = output.loc[output["0"] == "match",]
hold = hold.loc[hold["2"] == "M",]
for i in range(len(hold)):
    if hold.iat[i,1] == "K22955" or hold.iat[i,1] == "K13061" or hold.iat[i,1] == "K13061(2)" or hold.iat[i,1] == "K13060" or hold.iat[i,1] == "K22955(2)":
        acp = acp + 1
    elif hold.iat[i,1] == "K18096" or hold.iat[i,1] == "18096(2)":
        coa = coa + 1
    else:
        n = n + 1
total = acp + coa + n
mob["match"] = [acp, coa, n, total]
#unmatched stats
acp = 0
coa = 0
n = 0
hold = output.loc[output["0"] == "inconsistent",]
hold = hold.loc[hold["2"] == "M",]
for i in range(len(hold)):
    if hold.iat[i,1] == "K22955" or hold.iat[i,1] == "K13061" or hold.iat[i,1] == "K13061(2)" or hold.iat[i,1] == "K13060" or hold.iat[i,1] == "K22955(2)":
        acp = acp + 1
    elif hold.iat[i,1] == "K18096" or hold.iat[i,1] == "18096(2)":
        coa = coa + 1
    else:
        n = n + 1
total = acp + coa + n
mob["inconsistent"] = [acp, coa, n, total]
#ND stats
acp = 0
coa = 0
n = 0
hold = output.loc[output["0"] == "ND",]
hold = hold.loc[hold["2"] == "M",]
for i in range(len(hold)):
    if hold.iat[i,1] == "K22955" or hold.iat[i,1] == "K13061" or hold.iat[i,1] == "K13061(2)" or hold.iat[i,1] == "K13060" or hold.iat[i,1] == "K22955(2)":
        acp = acp + 1
    elif hold.iat[i,1] == "K18096" or hold.iat[i,1] == "18096(2)":
        coa = coa + 1
    else:
        n = n + 1
total = acp + coa + n
mob["ND"] = [acp, coa, n, total]
#total stats
mob["total"] = mob.sum(axis = 1)
mob = mob/2
