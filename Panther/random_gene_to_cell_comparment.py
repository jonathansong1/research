#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 18:07:41 2019

@author: jonathansong
"""
import pandas as pd 
import numpy as np
import math

# read mutation_map and only keep essential columns
mutation_map = pd.read_csv('/Users/jonathansong/Research_Folder/mutation_map.tsv', sep='\t')
# mutation_map = pd.read_csv('/Users/jonathansong/Research_Folder/mutation_map_exon.tsv', sep='\t')
# data_vars = mutation_map.columns.values.tolist()
to_keep = ['seqname', 'start', 'end', 'gene_id', 'gene_name']
mutation_map = mutation_map[to_keep]

# subset mutation_map to make it more easily accessible by chr number when looping
mutation_map['seqname'] = mutation_map['seqname'].astype(str)
mutation_map_1 = mutation_map.loc[mutation_map['seqname'] == '1']
mutation_map_2 = mutation_map.loc[mutation_map['seqname'] == '2']
mutation_map_3 = mutation_map.loc[mutation_map['seqname'] == '3']
mutation_map_4 = mutation_map.loc[mutation_map['seqname'] == '4']
mutation_map_5 = mutation_map.loc[mutation_map['seqname'] == '5']
mutation_map_6 = mutation_map.loc[mutation_map['seqname'] == '6']
mutation_map_7 = mutation_map.loc[mutation_map['seqname'] == '7']
mutation_map_8 = mutation_map.loc[mutation_map['seqname'] == '8']
mutation_map_9 = mutation_map.loc[mutation_map['seqname'] == '9']
mutation_map_10 = mutation_map.loc[mutation_map['seqname'] == '10']
mutation_map_11 = mutation_map.loc[mutation_map['seqname'] == '11']
mutation_map_12 = mutation_map.loc[mutation_map['seqname'] == '12']
mutation_map_13 = mutation_map.loc[mutation_map['seqname'] == '13']
mutation_map_14 = mutation_map.loc[mutation_map['seqname'] == '14']
mutation_map_15 = mutation_map.loc[mutation_map['seqname'] == '15']
mutation_map_16 = mutation_map.loc[mutation_map['seqname'] == '16']
mutation_map_17 = mutation_map.loc[mutation_map['seqname'] == '17']
mutation_map_18 = mutation_map.loc[mutation_map['seqname'] == '18']
mutation_map_19 = mutation_map.loc[mutation_map['seqname'] == '19']
mutation_map_20 = mutation_map.loc[mutation_map['seqname'] == '20']
mutation_map_21 = mutation_map.loc[mutation_map['seqname'] == '21']
mutation_map_22 = mutation_map.loc[mutation_map['seqname'] == '22']
mutation_map_X = mutation_map.loc[mutation_map['seqname'] == 'X']
mutation_map_Y = mutation_map.loc[mutation_map['seqname'] == 'Y']
mutation_map_list = [mutation_map_1, 
                     mutation_map_2,
                     mutation_map_3,
                     mutation_map_4,
                     mutation_map_5,
                     mutation_map_6,
                     mutation_map_7,
                     mutation_map_8,
                     mutation_map_9,
                     mutation_map_10,
                     mutation_map_11,
                     mutation_map_12,
                     mutation_map_13,
                     mutation_map_14,
                     mutation_map_15,
                     mutation_map_16,
                     mutation_map_17,
                     mutation_map_18,
                     mutation_map_19,
                     mutation_map_20,
                     mutation_map_21,
                     mutation_map_22,
                     mutation_map_X,
                     mutation_map_Y
                     ]
final_map_list = []

for mut_map in mutation_map_list:
    n = len(mut_map.index)
    bin_size = math.sqrt(n)
    n_bins = math.ceil(n/bin_size)
    final_map_list.append(np.array_split(mut_map, n_bins))

# construct comppi dictionary
# comppi_dict key:ensemble ID, value: locations    
file_path = "/Users/jonathansong/Research_Folder/Panther/comppi--proteins_locs--tax_hsapiens_loc_all.txt"
comppi = pd.read_csv(file_path, sep='\t')
key_list = []
value_list = []
comppi_dict = {}
for idx,row in comppi.iterrows():
    temp_list = row['Synonyms'].split('|')
    protein_name = row['Protein Name']
    for synonym in temp_list:
        if synonym[:4] == 'ENSG':
            if synonym not in comppi_dict:
                comppi_dict[synonym] = {protein_name: comppi.at[idx, 'Major Loc With Loc Score']}
            else:
                comppi_dict[synonym][protein_name] = comppi.at[idx, 'Major Loc With Loc Score']

# add new columns for incoming data
file_path = "/Users/jonathansong/Downloads/HUMAN_9606_idmapping.txt"
f = open(file_path, "r")
ID_list = []

for line in f:
    splitlist = line.split()
    if splitlist[1] == "Ensembl":
        ID_list.append(splitlist[2])
        
f.close()

#random_ID = []
#n = 200
#for i in range (0,n):
#    random_ID.append(random.choice(ID_list))

data = pd.DataFrame({'Ensemble_ID': ID_list})
major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
for item in major_loc_list:
    data[item] = "0"

    
# data = data[data['Validation_outcome'] == '1']

# bad_list contains genes not found in comppi
bad_chrompos_list = []

for idx,row in data.iterrows():
    try:
        protein_count = 0
        for loc_and_score_multiple in comppi_dict[data.at[idx, 'Ensemble_ID']].values():
            protein_count += 1
            temp_list = loc_and_score_multiple.split('|')
            for loc_and_score in temp_list:
                loc,score = loc_and_score.split(':')
                data.at[idx, loc] = str(float(data.at[idx, loc]) + float(score))
#                data_avg.at[idx, loc] = str(float(data_avg.at[idx, loc]) + float(score))
        for location in major_loc_list:
#            data_avg.at[idx, location] = str(float(data_avg.at[idx, location])/float(protein_count))
            data.at[idx, location] = str(float(data.at[idx, location])/float(protein_count))

    except:
        bad_chrompos_list.append(data.at[idx, 'Ensemble_ID'])
    if idx%2000 == 0:
        print(idx)
    
data.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/random_location_data.tsv')


