# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:05:35 2022

@author: mudoe
"""

#%%imports
import pandas as pd

#%%search for HSL
std = pd.read_csv("c:/users/mudoe/desktop/hetexp working/EGM HSL Discount.csv")
q = "jlw10r"
query = pd.read_csv("c:/users/mudoe/desktop/feature lists/" + q + ".cdf_ADAPchromatograms_deconvoluted1_adductless.txt.csv")
query2 = pd.read_csv("c:/users/mudoe/desktop/feature lists/" + q + "20.cdf_ADAPchromatograms_deconvoluted1_adductless.txt.csv")
hit = []
for i in range(len(std)):
    v = query.where(query[q+".cdf Peak m/z"] < std.iat[i,1]+.3)
    v = v.where(std.iat[i,1]-.3 < query[q+".cdf Peak m/z"])
    v = v.where(v["row retention time"] > std.iat[i,2]-.1)
    v = v.where(v["row retention time"] < std.iat[i,2]+.1)
    check = v.dropna(subset = ["row retention time"])
    if len(check) > 0:
        check = check.sort_values(by = [q+".cdf Peak height"], ascending=False)
        hit.append([std.iat[i,0],"{:.2e}".format(check.iat[0,2])])

#%%search for 102
std = pd.read_csv("c:/users/mudoe/desktop/hetexp working/EGM HSL Discount.csv")
#query = pd.read_csv("c:/users/mudoe/desktop/feature lists/" + q + ".cdf_ADAPchromatograms_deconvoluted1_adductless.txt.csv")
query2 = pd.read_csv("c:/users/mudoe/desktop/feature lists/" + q + "20.cdf_ADAPchromatograms_deconvoluted1_adductless.txt.csv")
hit2 = []
for i in range(len(std)):
    v = query2.where(query2[q+"20.cdf Peak m/z"] < 102.3)
    v = v.where(101.7 < query2[q+"20.cdf Peak m/z"])
    v = v.where(v["row retention time"] > std.iat[i,2]-.1)
    v = v.where(v["row retention time"] < std.iat[i,2]+.1)
    check2 = v.dropna(subset = ["row retention time"])
    if len(check2) > 0:
        check2 = check2.sort_values(by = [q+"20.cdf Peak height"], ascending=False)
        if check2.iat[0,2] > 25000:
            hit2.append([std.iat[i,0],"{:.2e}".format(check2.iat[0,2])])
