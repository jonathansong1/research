#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 13:53:22 2019

@author: jonathansong
"""
import pandas as pd
import numpy as np

raw_tsv_name = "/Users/jonathansong/Research_Folder/Immunodataset for Jonathan  - Sheet2.tsv"
data = pd.read_csv(raw_tsv_name, sep='\t')

validation_column = "Validation_outcome"

print(data[validation_column].value_counts())

to_keep = ['Identifier',
           'Patient_ID',
        'Gene_symbol',
         'Mutant_peptide',
         'Wildtype_peptide',
         'AA_mutation',
         'Chrom_pos',
         'Validation_outcome']

data = data[to_keep]

data.dropna(axis=0, inplace=True, how='any')

data[validation_column]=np.where(data[validation_column] == 'YES', "1", data[validation_column])
data[validation_column]=np.where(data[validation_column] == 'NO', "0", data[validation_column])
data[validation_column]=np.where(data[validation_column] == 'Y', "1", data[validation_column])
data[validation_column]=np.where(data[validation_column] == 'N', "0", data[validation_column])

yes_no_drop_list = []
for idx, row in data.iterrows():
    if data.at[idx, validation_column] != "1" and data.at[idx, validation_column] != "0":
        yes_no_drop_list.append(idx)
data.drop(labels=yes_no_drop_list, axis=0, inplace=True)

print(data[validation_column].value_counts())

"""
# to make dataset have equal #'s of 0 and 1
import random
negative_bias_drop_list = []
counts = data[validation_column].value_counts()
cutoff_ratio = counts[1] / counts[0]
for idx, row in data.iterrows():
    if data.at[idx, validation_column] == "0":
        if random.random() > cutoff_ratio:
            negative_bias_drop_list.append(idx)
data.drop(labels=negative_bias_drop_list, axis=0, inplace=True)
print(data[validation_column].value_counts())
"""
# data["Panther_gene_symbol"] = data["Gene_symbol"]

data.to_csv(sep='\t', path_or_buf='clean_data.tsv')
