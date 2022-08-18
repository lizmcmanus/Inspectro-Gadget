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
    GABAA_Alpha = np.nansum(data['GABAA_Alpha'])
    GABAA_Beta = np.nansum(data['GABAA_Beta'])
    GABAA_Gamma = np.nansum(data['GABAA_Gamma'])
    
    GABAA = GABAA_Alpha + GABAA_Beta + GABAA_Gamma
    
    #calculating excitation-inhibiiton ratio
    EI_ratio = (AMPA+NMDA)/GABAA
    return EI_ratio
     