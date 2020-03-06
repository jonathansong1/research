#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:57:48 2019

@author: jonathansong
"""
import pandas as pd

data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")    

mutant_peptide_column_name = "Mutant_peptide"

tiny_list = ["A", "C", "G", "S", "T"]
small_list = ["A", "B", "C", "D", "G", "N", "P", "S", "T", "V"]
aliphatic_list = ["A", "I", "L", "V"]
aromatic_list = ["F", "H", "W", "Y"]
nonpolar_list = ["A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y"]
polar_list = ["D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z"]
charged_list = ["B", "D", "E", "H", "K", "R", "Z"]
basic_list = ["H", "K", "R"]
acidic_list = ["B", "D", "E", "Z"]
sulfur_list = ["C", "M"]

list_list = [tiny_list, small_list, aliphatic_list, aromatic_list,
             nonpolar_list, polar_list, charged_list, basic_list,
             acidic_list, sulfur_list]

data["Tiny.mol.p"] = ""
data["Small.mol.p"] = ""
data["Aliphatic.mol.p"] = ""
data["Aromatic.mol.p"] = ""
data["Nonpolar.mol.p"] = ""
data["Polar.mol.p"] = ""
data["Charged.mol.p"] = ""
data["Basic.mol.p"] = ""
data["Acidic.mol.p"] = ""
data["Sulfur.mol.p"] = ""

column_name_list = ["Tiny.mol.p", "Small.mol.p", "Aliphatic.mol.p",
                    "Aromatic.mol.p", "Nonpolar.mol.p", "Polar.mol.p",
                    "Charged.mol.p", "Basic.mol.p", "Acidic.mol.p",
                    "Sulfur.mol.p"]

for idx, row in data.iterrows():
    divide = len(data.at[idx, mutant_peptide_column_name])
    tiny_count = 0
    small_count = 0
    aliphatic_count = 0
    aromatic_count = 0
    nonpolar_count = 0
    polar_count = 0
    charged_count = 0
    basic_count = 0
    acidic_count = 0
    sulfur_count = 0
    count_list = [tiny_count, small_count, aliphatic_count, aromatic_count,
             nonpolar_count, polar_count, charged_count, basic_count,
             acidic_count, sulfur_count]
    for i in range(len(list_list)):
        for letter in list_list[i]:
            count_list[i] += data.at[idx, mutant_peptide_column_name].count(letter)
    for j in range(len(column_name_list)):
        data.at[idx, column_name_list[j]] = (count_list[j] / divide)


data.to_csv(sep='\t', path_or_buf='clean_data.tsv')