#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 13:13:32 2020

@author: jonathansong
"""

import pandas as pd

rizvi = pd.read_excel("/Users/jonathansong/Research_Folder/aaa1348_TableS5.xlsx", header=0)
data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")

data_subset = data[data["Identifier"] == "rizvi_PMID25765070"]

for idx,row in data_subset.iterrows():
    rizvi_subset = rizvi[rizvi["Gene"] == data_subset.at[idx, "Gene_symbol"]]
    for idx1,row1 in rizvi_subset.iterrows():
        begin = rizvi_subset.at[idx1, "Ref"]
        end = rizvi_subset.at[idx1, "Alt"]
        if type(begin) == str and type(end) == str:
            data.at[idx, "allele"] = begin + "/" + end
            data_subset.at[idx, "allele"] = begin + "/" + end
            
data.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/clean_data.tsv')
