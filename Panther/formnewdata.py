#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 13:49:04 2019

@author: jonathansong
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 13:53:22 2019

@author: jonathansong
"""
import pandas as pd
import numpy as np

raw_tsv_name = "Immunodataset for Jonathan  - Sheet2 copy.tsv"
data = pd.read_csv(raw_tsv_name, sep='\t')

validation_column = "Validation_outcome"

print(data[validation_column].value_counts())

to_keep = ['Gene_symbol',
         'Mutant_peptide',
         'Validation_outcome',
         'Panther_gene_symbol']

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

data_1 = data[data[validation_column] == "1"]
data_0 = data[data[validation_column] == "0"]

data_1.to_csv(sep='\t', path_or_buf='1_panther_test_data.tsv')
data_0.to_csv(sep='\t', path_or_buf='0_panther_test_data.tsv')

file = open("1_to_panther.txt","w+")
for idx, row in data_1.iterrows():
    file.write(data_1.at[idx, "Gene_symbol"])
    file.write("\n")
file.close()

file = open("0_to_panther.txt","w+")
for idx, row in data_0.iterrows():
    file.write(data_0.at[idx, "Gene_symbol"])
    file.write("\n")
file.close()

file = open("1and0_to_panther.txt","w+")
symbol_list = []
for idx_1, row_1 in data_1.iterrows():
    for idx_0, row_0 in data_0.iterrows():
        if data_1.at[idx_1, "Gene_symbol"] == data_0.at[idx_0, "Gene_symbol"]:
            file.write(data_1.at[idx_1, "Gene_symbol"])
            file.write("\n")
file.close()