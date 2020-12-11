# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Tool to show the receptor expressions within an MRS region

MRS regions masks must be (in MNI space?))

GABA and Glutamate regions must ine in a .tsv file "receptors.tsv"
"""
import os
import numpy as np
import pandas as pd
import nibabel as ni
from extract_mrna import extract_mrna
from make_violin import make_violin

# Project directory
data_dir = 'E://Taiwan/Inspectro-Gadget/'

# Directory for outputs
output_dir = os.path.join(data_dir,'expression')

# Create outputs directory if missing
#if not os.path.isdir(output_dir):
#  os.mkdir(output_dir)

# Load MRS region
region_mask = ni.load(os.path.join(data_dir,'MRS_region.nii.gz')).get_fdata()
mask_size = region_mask[region_mask==1]
mask_size = len(mask_size)

# Load in the list of receptors and subunits
receptors_list = pd.read_csv(data_dir+'receptors.tsv', delimiter='\t', header=None)
#extract unique receptors
receptors = receptors_list[1].unique()

# loop thorugh each receptor type
for x, receptor_t in enumerate(receptors):
    
    #makes variables for subunits of each specific receptor type
    receptor_type = receptors_list.loc[receptors_list[1] == receptor_t]
    receptor_region_data = np.zeros([mask_size,len(receptor_type[0])])
     
    #loop through subunits within each receptor type
    for a, sub_a in enumerate(receptor_type[0]):
        #call function to extract mrna data within each region
        extract_mrna(a,sub_a,receptor_region_data, region_mask)
    
    #function to make violin plots and save as PDF's
    make_violin(receptor_region_data , receptor_type, receptor_t)
  