#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 12:41:06 2020

@author: jonathansong
"""
import pandas as pd
data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")
f = open("/Users/jonathansong/Research_Folder/transvar.txt", "w+")


for idx,row in data.iterrows():
    if type(data.at[idx, "allele"]) != str:
        string = str(row["Gene_symbol"]) + ":p." + str(row["AA_mutation"])
        f.write(string + '\n')
    
f.close()

