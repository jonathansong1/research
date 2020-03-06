#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:37:45 2019

@author: jonathansong
"""
import pandas as pd

    
file_path = "/Users/jonathansong/Downloads/HUMAN_9606_idmapping.txt"
f = open(file_path, "r")
AC_list = []
GN_list = []

for line in f:
    splitlist = line.split()
#    AC_list.append(splitlist[0])
#    GN_list.append(splitlist[2])

    if splitlist[1] == "Gene_Name":
        AC_list.append(splitlist[0])
        GN_list.append(splitlist[2])
#    elif splitlist[1] == "Gene_Synonym":
#        AC_list.append(splitlist[0])
#        GN_list.append(splitlist[2])

#lists = {"Protein Name":AC_list, "Gene Name":GN_list}
#uniprot = pd.DataFrame(lists)

uniprot = dict(zip(GN_list, AC_list))


data = pd.read_csv("/Users/jonathansong/Research_Folder/out_data.tsv", sep="\t")
major_loc_list = ['N/A', 'membrane', 'cytosol', 'nucleus', 'secretory-pathway', 'mitochondrion', 'extracellular']
for item in major_loc_list:
    data[item] = "0"

file_path = "/Users/jonathansong/Research_Folder/Panther/comppi--proteins_locs--tax_hsapiens_loc_all.txt"
comppi = pd.read_csv(file_path, sep='\t')
#uniprot_id_list = []
#location_list = []
#for idx, row in comppi.iterrows():
#    uniprot_id_list.append(comppi.at[idx, 'Protein Name'])
#    location_list.append(comppi.at[idx, 'Major Loc With Loc Score'])
#    syn_list = comppi.at[idx, 'Synonyms']
#    for synonym in syn_list:
#        uniprot_id_list.append(synonym)
#        location_list.append(comppi.at[idx, 'Major Loc With Loc Score'])
#

    


bad_genes = []
multiple_genes = []

for idx, row in data.iterrows():
    try:
        gene_symbol = data.at[idx, 'Gene_symbol']
        protein_name = uniprot[gene_symbol]
        dict_subset = comppi.loc[comppi['Protein Name'] == protein_name]
        temp_list = (dict_subset['Major Loc With Loc Score'].iloc[0]).split('|')
        for loc in temp_list:
            loc = loc.split(':')
            data.at[idx, loc[0]] = loc[1]
    except:
        
        if gene_symbol in uniprot.keys():
            multiple_genes.append(gene_symbol)
        else:
            bad_genes.append(gene_symbol)
data.to_csv(sep='\t', path_or_buf='add_location_data.tsv')