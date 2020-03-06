#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:46:32 2019

@author: jonathansong
"""
import pandas as pd
import numpy as np


def N_and_C_terminus(data_input):
    
    data=data_input
    
    data["N.terminus"] = ""
    data["C.terminus"] = ""
    
    for idx, row in data.iterrows():
        data.at[idx, "N.terminus"] = data.at[idx, "Mutant.peptide"][0]
        data.at[idx, "C.terminus"] = data.at[idx, "Mutant.peptide"][-1]
        
    subset = pd.get_dummies(data["N.terminus"], prefix="N.terminus")
    subset = subset.join(pd.get_dummies(data["C.terminus"], prefix="C.terminus"))
    
    noise = subset + 0.001 * np.random.rand(subset.shape[0], subset.shape[1])
    
    del_vars = ["N.terminus",
                "C.terminus"]
    data_vars = data.columns.values.tolist()
    data_final_vars = [i for i in data_vars if i not in del_vars]
    data_final = data[data_final_vars]
    
    data_final_and_noise = data_final.join(noise)
    
    return data_final_and_noise
    

    