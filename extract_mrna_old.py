# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 10:25:43 2022

@author: lizmc
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:13:59 2020

@author: lizmc
"""
import os
import nibabel as ni
import numpy as np
from scipy.stats import iqr
# data_dir = 'E://Taiwan/Inspectro-Gadget/'

def extract_mrna(subunit_path, region_mask):
    mrna = ni.load(subunit_path).get_fdata()
    mrna_removed = mrna[(mrna != 0)]
    #finding max value for whole brain to normalise data once region has been extracted
    mrna_max = mrna_removed.max()
    mrna_mean = mrna_removed.mean()
    mrna_median = np.median(mrna_removed)
    #Indexing each expression by the mask- shows only the receptor expression within the MRS ROI
    region = mrna[region_mask==1]
    #normalising region mrna values by whole brain receptors
    region_norm = region/mrna_max 
    #region_norm = region/mrna_mean
    #array for all normalised values for all subunits

    return region_norm
