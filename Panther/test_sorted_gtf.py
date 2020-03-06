#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:12:03 2020

@author: jonathansong
"""

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
#mutation_map = pd.read_csv('/Users/jonathansong/Research_Folder/mutation_map.tsv', sep='\t')
# mutation_map = pd.read_csv('/Users/jonathansong/Research_Folder/mutation_map_exon.tsv', sep='\t')
mutation_map = pd.read_csv('/Users/jonathansong/Research_Folder/sorted_mutation_map.tsv', sep='\t')

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
data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")
major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
for item in major_loc_list:
    data[item] = '0'
new_col_list = ['Ensemble_ID', 'New_Gene_ID']
for item in new_col_list:
    data[item] = '0'
    
data['Validation_outcome'] = data['Validation_outcome'].astype(str)    

# bad_list contains genes not found in comppi
bad_chrompos_list = []
bad_sort_list = []
time_list = []
import time
## set starting time
#start_time = time.time()
## print time since start to finding ensemble ID
#print("--- %s seconds ---" % (time.time() - start_time))

for idx,row in data.iterrows():
    start_time = time.time()
    chrom, mut_loc = row['Chrom_pos'].split(':')
    chrom = chrom[3:]
    if chrom == 'X':
        chrom = '23'
    elif chrom == 'Y':
        chrom = '24'
        
    # select correct subset of mutation_map
    
    mutation_map_chrom = final_map_list[int(chrom)-1]
    sortworked = False
    for i in range(0, len(mutation_map_chrom)):
        map_n = mutation_map_chrom[i]
        map_len = len(map_n.index)
        if int(mut_loc) >= int(map_n.iloc[0]['start']):
            if int(mut_loc) <= int(map_n.iloc[map_len - 1]['end']):
                working_map = map_n
                sortworked = True
                break
    if not sortworked:
        bad_sort_list.append(idx)
    
    # loop through mutation_map subset to find where mut_loc is in the range of some gene
    # save the ensemble id and gene name
    for idx_mut, row_mut in working_map.iterrows():
        if int(mut_loc) >= int(row_mut['start']):
            if int(mut_loc) <= int(row_mut['end']):
                data.at[idx, 'Ensemble_ID'] = row_mut['gene_id']
                data.at[idx, 'New_Gene_ID'] = row_mut['gene_name']
                break
    try:
        protein_count = 0;
        loc_and_score_list = comppi_dict[data.at[idx, 'Ensemble_ID']].values()
        for multiple_loc_and_score in loc_and_score_list:
            protein_count += 1
            temp_list = multiple_loc_and_score.split('|')
            for loc_and_score in temp_list:
                loc,score = loc_and_score.split(':')
                data.at[idx, loc] = str(float(data.at[idx, loc]) + float(score))
        for location in major_loc_list:
            data.at[idx, location] = str(float(data.at[idx, location])/float(protein_count))
    except:
        bad_chrompos_list.append(idx)
    
    time_list.append(time.time() - start_time)    
    if idx%10 == 0:
        print(idx)

total_time = sum(time_list)
print('Total Run Time: ', total_time) 
print('Average Run Time Per Peptide', total_time/len(time_list))      
# data.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/clean_data.tsv')
# data_avg.to_csv(sep='\t', path_or_buf='add_avg_location_data.tsv')

not_there = []
for index in bad_chrompos_list:
    chrom, mut_loc = data.at[index, 'Chrom_pos'].split(':')
    chrom = chrom[3:]
    if chrom == 'X':
        chrom = '23'
    elif chrom == 'Y':
        chrom = '24'
#    print(chrom)
#    print(mut_loc)
mutation_map_chrom = final_map_list[int(chrom)-1]
sortworked = False
for i in range(0, len(mutation_map_chrom)):
    working_map = mutation_map_chrom[i]
    for idx_mut, row_mut in working_map.iterrows():
        if int(mut_loc) >= int(row_mut['start']):
            if int(mut_loc) <= int(row_mut['end']):
                if data.at[idx, 'Ensemble_ID'] != row_mut['gene_id']:
                    print(row_mut['gene_id'])
                    print(row_mut['gene_name'])
#                    print(row_mut['start'])
#                    print(row_mut['end'])
#                    print(i)
#                    not_there.append(idx)