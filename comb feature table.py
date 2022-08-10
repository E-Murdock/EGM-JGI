# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 11:52:14 2022

@author: mudoe
"""
#may want to update std2.csv as the other ones have changed a bit

#%% imports
import os as os
import pandas as pd
import numpy as np
import math
import copy
os.chdir("c:/users/mudoe/desktop/hetexp working")

#pull data
std = pd.read_csv("std2.csv")
#these standards currently neglect C4, 3OHC4
qu = pd.read_csv("Aaron Query.csv")
ftfull = pd.read_csv("full.ft.csv")
alp = "20211117_JGI-AK-TH_AP_507288_EColiHSLs_final_QE-139_C18_USDAY63649_POS_MSMS_369_EAWP156-NEG-3OHC10HSL_A_Rg80to1200-CE102040-MediaNegCtrl-Bl21-S1_Run1152.mzML Peak height"
bet = "20211117_JGI-AK-TH_AP_507288_EColiHSLs_final_QE-139_C18_USDAY63649_POS_MSMS_370_EAWP156-NEG-3OHC10HSL_B_Rg80to1200-CE205060-MediaNegCtrl-Bl21-S1_Run1159.mzML Peak height"
gam = "20211117_JGI-AK-TH_AP_507288_EColiHSLs_final_QE-139_C18_USDAY63649_POS_MSMS_371_Media-control_A_Rg80to1200-CE102040-MediaCtrl-NoInoc-S1_Run1149.mzML Peak height"
delt = "20211117_JGI-AK-TH_AP_507288_EColiHSLs_final_QE-139_C18_USDAY63649_POS_MSMS_372_Media-control_B_Rg80to1200-CE205060-MediaCtrl-NoInoc-S1_Run1162.mzML Peak height"
ft = ftfull[ftfull['row ID'].isin(qu["_dyn_#scan"])]
Na = 22.989218
H = 1.007276

#%% finding M+H
o = []
p = []
#for i in range(18):
for i in range(len(std)):
    v = ft.where(ft["row m/z"] > std.iat[i,3]-(std.iat[i,3]/200000))
    v = v.where(v["row m/z"] < std.iat[i,3]+(std.iat[i,3]/200000))
    v = v.where(v["row retention time"] > std.iat[i,2]-.1)
    v = v.where(v["row retention time"] < std.iat[i,2]+.1)
    check = v.dropna(subset = ["row ID"])
#this is the section you'll want to change if you go back to averages  
    if len(check) == 1:
        temp = np.average(check, axis=0)
        o.append(temp)
        p.append(std.iat[i,0])
    if len(check) > 1:
        q = []
        check = check.reset_index()
        for j in range(len(check)):
            q.append(abs(std.iat[i,3]-check.iat[j,3]))
            z = min(q)
            z = q.index(z)
        temp = pd.DataFrame(columns=ft.columns)
        temp.loc[0,:] = check.loc[z,:]
        temp = np.average(temp, axis=0)
        o.append(temp)
        p.append(std.iat[i,0])
##see above
hknowns = pd.DataFrame(index=range(len(o)),columns=ft.columns)
for i in range(len(o)):
    hknowns.loc[i,:] = o[i]
hknowns["row ID"] = p
#hknowns = hknowns.drop("Unnamed: 415",axis=1)

#%% finding M+Na
o = []
p = []
for i in range(len(std)):
    v = ftfull.where(ftfull["row m/z"] > (std.iat[i,3]-H+Na)-((std.iat[i,3]-H+Na)/100000))
    v = v.where(v["row m/z"] < (std.iat[i,3]-H+Na)+((std.iat[i,3]-H+Na)/100000))
    v = v.where(v["row retention time"] > std.iat[i,2]-.11)
    v = v.where(v["row retention time"] < std.iat[i,2]+.11)
    check = v.dropna(subset = ["row ID"])
#this is the section you'll want to change if you go back to averages  
    if len(check) == 1:
        temp = np.average(check, axis=0)
        o.append(temp)
        p.append(std.iat[i,0])
    if len(check) > 1:
        q = []
        check = check.reset_index()
        for j in range(len(check)):
            q.append(abs(hknowns.iat[i,3]-check.iat[j,3]))
        z = min(q)
        z = q.index(z)
        temp = pd.DataFrame(columns=ft.columns)
        temp.loc[0,:] = check.loc[z,:]
        temp = np.average(temp, axis=0)
        o.append(temp)
        p.append(std.iat[i,0])
##see above
naknowns = pd.DataFrame(index=range(len(o)),columns=ft.columns)
for i in range(len(o)):
    naknowns.loc[i,:] = o[i]
naknowns["row ID"] = p
#naknowns = naknowns.drop("Unnamed: 415",axis=1)
naknowns = naknowns[naknowns["row ID"].isin(hknowns["row ID"])]

#%% verify the signals within the knowns
background_ratio = 3
control = hknowns[[alp,bet,gam,delt]]
control = control.max(axis=1)
for j in range(3,len(hknowns.columns)):
    for i in range(len(hknowns)):
        if hknowns.iat[i,j] < (background_ratio*control.iat[i]):
            hknowns.iat[i,j] = 0
control = naknowns[[alp,bet,gam,delt]]
control = control.max(axis=1)
for j in range(3,len(naknowns.columns)):
    for i in range(len(naknowns)):
        if naknowns.iat[i,j] < (background_ratio*control.iat[i]):
            naknowns.iat[i,j] = 0
for i in range(len(hknowns)):
    for j in range(3,hknowns.shape[1],1):
        if (.1*hknowns.iat[i,j]) > naknowns.iat[i,j]:
        #if hknowns.iat[i,j] > ((naknowns.iat[i,j])-((naknowns.iat[i,j])/10)):
            hknowns.iat[i,j] = 0
        if hknowns.iat[i,j] < (10**6):
            hknowns.iat[i,j] = 0
        if naknowns.iat[i,j] < (10**6):
            hknowns.iat[i,j] = 0
s = copy.deepcopy(hknowns)
            
#%% order runs and check for spillover
hknowns2 = copy.deepcopy(hknowns)
cols = hknowns2.columns
colsrun=["row ID","row m/z","row retention time"]
#cols2 = []
for i in range(3,len(cols)):
    x = cols[i]
    x = x.split("Run",1)[1]
    x = x.split(".",1)[0]
    colsrun.append(int(x))
hknowns2.columns = colsrun
cols3=colsrun[3:415]
cols3.sort()
colsrun[3:415]=cols3
t = hknowns2[["row ID","row m/z","row retention time"]]
hknowns2 = hknowns2.drop("row ID",axis=1)
hknowns2 = hknowns2.drop("row m/z",axis=1)
hknowns2 = hknowns2.drop("row retention time",axis=1)
hknowns2 = hknowns2.sort_index(axis=1,level=None,ascending=True,inplace=False)
t = pd.concat([t, hknowns2], axis=1)
q = []
q2=[]
for j in range(4,len(colsrun)): #note 4 not 3 because 8 is first run
    for i in range(len(t)):
        if t.iat[i,j] > (10**5):
            if t.iat[i,j] < (10**7):
                if t.iat[i,j-1] > (10**8):
                    if t.iat[i,j] < (t.iat[i,j-1]):
                        q2.append((cols3[j],cols3[j-1],t.iat[i,0],math.log10(t.iat[i,j]),math.log10(t.iat[i,j-1])))
                        q.append((colsrun[j],colsrun[j-1],t.iat[i,0],t.iat[i,j],t.iat[i,j-1]))
    
#%% get the locus tag, run number, and organism name
###
###t is runs, u is locus tags
###

colsloc = ["row ID","row m/z","row retention time"]
for i in range(3,len(cols)):
    x = cols[i]
    if x.find("Batch304-") == -1:
        if x.find("Std") == -1:
            colsloc.append(x)
    if x.find("Batch304-") != -1:
        x = x.split("Batch304-",1)[1]
        x = x.split("_",1)[0]
        colsloc.append(x)
    if x.find("Std") != -1:
        x = x.split("Std",1)[1]
        x = x.split("_",1)[1]
        x = x.split("_",1)[0]
        colsloc.append(x+" HSL standard")

user = pd.read_csv("user.csv")
colsloc4=copy.deepcopy(colsloc)
for i in range(len(user)):
    user.iat[i,0] = user.iat[i,0].split("_",1)[1]
for i in range(len(colsloc)):
    for j in range(len(user)):
        if colsloc[i] == user.iat[j,0]:
            colsloc[i] = user.iat[j,1]
locus = pd.read_csv("locus.csv")
colsloc = pd.Series(colsloc).to_frame()
for i in range(3,len(colsloc)):
    for j in range(len(locus)):
        if colsloc.iat[i,0] == locus.iat[j,2]:
            colsloc.iat[i,0] = locus.iat[j,0]
u = copy.deepcopy(hknowns)
u.columns = colsloc[0]
