#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 22:27:55 2020

@author: jonathansong
"""
import pandas as pd
import numpy as np

#file_path = "/Users/jonathansong/Research_Folder/Panther/comppi--proteins_locs--tax_hsapiens_loc_all.txt"
#comppi = pd.read_csv(file_path, sep='\t')
#df = pd.DataFrame()
#df["Protein Name"] = ""
#df["Ensemble ID"] = ""
#major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
#for loc in major_loc_list:
#    df[loc] = 0
#for idx,row in comppi.iterrows():
#    temp_list = row['Synonyms'].split('|')
#    for synonym in temp_list:
#        if synonym[:4] == 'ENSG':
#            temp_df = pd.DataFrame({"Protein Name" : [synonym]})
#            temp_df["Protein Name"] = row['Protein Name']
#            temp_df["Ensemble ID"] = synonym
#            major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
#            for loc in major_loc_list:
#                temp_df[loc] = 0
#            loc_list = row["Major Loc With Loc Score"].split('|')
#            for loc_and_score in loc_list:
#                loc,score = loc_and_score.split(':')
#                temp_df[loc] = float(score)
#            df = df.append(temp_df)
#        
#        
#
#f = open(file_path, "r")
#complete_reviewed_list = []
#i = 0
#for line in f:
#    temp = line.split()
#    if temp[0] == 'AC':
#        i += 1
#        temp.pop(0)
#        if i%2000 == 0:
#            print(i)
#        for item in temp:
#            item = item.strip(";")
#            complete_reviewed_list.append(item)
#f.close()

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

comppi_subset = comppi[comppi["Naming Convention"] == "UniProtKB/Swiss-Prot/P"]
complete_reviewed_list_2 = []
for idx,row in comppi_subset.iterrows():
    temp_list = row["Synonyms"].split("|")
    complete_reviewed_list_2.append(temp_list[0])

df = pd.DataFrame()
df["Ensemble ID"] = ""
df["Gene_symbol"] = ""
df["Protein Name"] = ""

major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
for loc in major_loc_list:
    df[loc] = 0

data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")

for idx,row in data.iterrows():
    if data.at[idx, 'Ensemble_ID'] in comppi_dict.keys():
        reviewed_list = []
        unreviewed_list = []
        reviewed = False
        for gene_id in comppi_dict[data.at[idx, 'Ensemble_ID']].keys():
            if gene_id in complete_reviewed_list_2:
                reviewed = True
                reviewed_list.append(gene_id)
            else:
                unreviewed_list.append(gene_id)
        if reviewed:
            gene_list = reviewed_list
        else:
            gene_list = unreviewed_list
        
        working_dict = comppi_dict[data.at[idx, 'Ensemble_ID']]
        for gene in gene_list:
            temp_df = pd.DataFrame({"Gene_symbol" : [data.at[idx, 'Gene_symbol']]})
            temp_df["Ensemble ID"] = data.at[idx, 'Ensemble_ID']
            temp_df["Protein Name"] = gene
            major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
            for loc in major_loc_list:
                temp_df[loc] = 0
            temp_list = working_dict[gene].split('|')
            for loc_and_score in temp_list:
                loc,score = loc_and_score.split(':')
                temp_df[loc] = score
            df = df.append(temp_df)
            
df.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/reviewed_unreviewed.tsv')

df = pd.read_csv('/Users/jonathansong/Research_Folder/reviewed_unreviewed.tsv', sep='\t')

for idx,row in df.iterrows():
    if df.at[idx, "Protein Name"] in complete_reviewed_list_2:
        df.at[idx, "Reviewed"] = "Reviewed"
    else:
        df.at[idx, "Reviewed"] = "Unreviewed"

df.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/reviewed_unreviewed.tsv')
df = pd.read_csv('/Users/jonathansong/Research_Folder/reviewed_unreviewed.tsv', sep='\t')

unique_keys, indices = np.unique(df["Protein Name"], return_index = True)
df = df[indices]
