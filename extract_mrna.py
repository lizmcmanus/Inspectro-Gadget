# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:13:59 2020

@author: lizmc
"""
import os
import nibabel as ni
import numpy as np
import math
# data_dir = 'E://Taiwan/Inspectro-Gadget/'

def extract_mrna(subunit_path, region_mask):
    mrna = ni.load(subunit_path).get_fdata()
    mrna_removed = np.array(mrna[(mrna != 0)])
    

    #robust sigmoid IQR
    median = np.median(mrna_removed)
    Q1 = np.percentile(mrna_removed, 25)
    Q3 = np.percentile(mrna_removed, 75)
    IQRx = (Q3-Q1)/1.35
    
    
    #mrna_removed = mrna[(mrna != 0)]
    
    #robust sigmoid
    Xy = 1/(1+ np.exp(-(mrna-median)/IQRx))
    region_norm = Xy[region_mask==1]
    
    
   
    

    
    

   

    return region_norm