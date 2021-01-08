# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 13:14:26 2021

@author: lizmc
"""
import operator
import numpy as np 
def ex_in(data):
    #summing array for each subunit
    AMPA = np.nansum(data['AMPA'])
    NMDA = np.nansum(data['NMDA'])
    GABAa = np.nansum(data['GABAA'])
    
    #calculating excitation-inhibiiton ratio
    EI_ratio = (AMPA+NMDA)/GABAa
    return EI_ratio
     