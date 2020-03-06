#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 12:54:18 2020

@author: jonathansong
"""

#import pandas as pd
#
#data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")
#
#
#bad_seq = []
#data["correction_AA_mutation"] = ""
#for idx,row in data.iterrows():
#    mut = data.at[idx, "Mutant_peptide"]
#    wild = data.at[idx, "Wildtype_peptide"]
#    if len(mut) == len(wild):
#        for i in range(0, len(mut)):
#            if mut[i] != wild[i]:
#                loc = data.at[idx, "AA_mutation"][1:-1]
#                output = str(wild[i]) + loc + str(mut[i])
#                if output != data.at[idx, "AA_mutation"]:
#                    data.at[idx, "correction_AA_mutation"] = output
#                
#    else:
#        bad_seq.append(idx)

import pandas as pd

data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")

bad_seq = []
for idx,row in data.iterrows():
    mut = data.at[idx, "Mutant_peptide"]
    wild = data.at[idx, "Wildtype_peptide"]
    if len(mut) == len(wild):
        for i in range(0, len(mut)):
            if mut[i] != wild[i]:
                loc = data.at[idx, "AA_mutation"][1:-1]
                output = str(wild[i]) + loc + str(mut[i])
                if output != data.at[idx, "AA_mutation"]:
                    data.at[idx, "AA_mutation"] = output
                
    else:
        bad_seq.append(idx)
    