# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 12:41:04 2022

@author: lizmc
"""

import operator
import numpy as np 

def ex_in_multi(data):
    #making an array for each recetor type
    AMPA = data[["GRIA1", "GRIA2", "GRIA3",	"GRIA4"]]
    NMDA = data[["GRIN1", "GRIN2A", "GRIN2B", "GRIN2C"]]
    GABAA = data[["GABRA1",	"GABRA2",	"GABRA3", "GABRA4",	"GABRA5",	"GABRB1",	"GABRB2",	"GABRB3",	"GABRG1",	"GABRG2", "GABRG3"]]
    
    #summing array for each subunit
    AMPA = AMPA.sum(axis=1) #sum each column
    AMPA = AMPA.sum(axis=0) #sum all together
    NMDA = NMDA.sum(axis=1) #sum each column
    NMDA = NMDA.sum(axis=0) #sum all together
    GABAA = GABAA.sum(axis=1) #sum each column
    GABAA = GABAA.sum(axis=0) #sum all together
    
    
    #calculating excitation-inhibiiton ratio
    EI_ratio = (AMPA+NMDA)/GABAA
    return EI_ratio