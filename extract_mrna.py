# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:13:59 2020

@author: lizmc
"""
import os
import glob
import nibabel as ni
data_dir = 'E://Taiwan/Inspectro-Gadget/'

def extract_mrna(subunit, region_mask):
    subunit_path = os.path.join(data_dir,'mRNA_images',subunit+'_mirr_mRNA.nii.gz')
    mrna = ni.load(subunit_path).get_fdata()
    #finding max value for whole brain to normalise data once region has been extracted
    mrna_max = mrna.max()
    #Indexing each expression by the mask- shows only the receptor expression within the MRS ROI
    region = mrna[region_mask==1]
    #normalising region mrna values by whole brain receptors
    region_norm = region/mrna_max
    #array for all normalised values for all GABAa subunits
    return region_norm
