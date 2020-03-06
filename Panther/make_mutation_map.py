#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 19:24:18 2019

@author: jonathansong
"""
from gtfparse import read_gtf
import pandas as pd

mutation_map = read_gtf("/Users/jonathansong/Downloads/Homo_sapiens.GRCh37.75.gtf")
mutation_map = mutation_map.loc[mutation_map['feature'] == 'exon']


to_keep = ['seqname', 'start', 'end', 'gene_id', 'gene_name']
mutation_map = mutation_map[to_keep]

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
print("done separating maps")
count = 0
for i in range(0, len(mutation_map_list)):
    mutation_map_list[i].sort_values(by=['start', 'end'], inplace=True)
    print("map", i+1, "done")
          
final_mutation_map = pd.concat(mutation_map_list)
print("concatenation done")
final_mutation_map.to_csv(sep='\t', path_or_buf='/Users/jonathansong/Research_Folder/mutation_map_exon.tsv')

    