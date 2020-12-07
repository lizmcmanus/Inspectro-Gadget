# -*- coding: utf-8 -*-
"""
Tool to show the receptor expressions within an MRS region

MRS regions masks must be (in MNI space?))

GABA and Glutamate regions must ine in a .tsv file "receptors.tsv"
"""

import os
import numpy as np
import pandas as pd
import glob
import nibabel as ni
import seaborn as sns

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

# Load in the list of subunits for each receptor types
GABAa = pd.read_csv(data_dir+'GABAa.tsv', delimiter='\t', header=None)[0]
n_GABAa = len(GABAa)

GABAb = pd.read_csv(data_dir+'GABAb.tsv', delimiter='\t', header=None)[0]
n_GABAb = len(GABAb)

AMPA = pd.read_csv(data_dir+'AMPA.tsv', delimiter='\t', header=None)[0]
n_AMPA = len(AMPA)

NMDA = pd.read_csv(data_dir+'NMDA.tsv', delimiter='\t', header=None)[0]
n_NMDA = len(NMDA)

mGlu = pd.read_csv(data_dir+'mGlu.tsv', delimiter='\t', header=None)[0]
n_mGlu = len(mGlu)

kainate = pd.read_csv(data_dir+'kainate.tsv', delimiter='\t', header=None)[0]
n_kainate = len(kainate)

# Load data for each GABAa subunit
GABAa_region_data = np.zeros([mask_size,n_GABAa])
for a, sub_a in enumerate(GABAa):
    GABAa_path = glob.glob(os.path.join(data_dir,'mRNA_Expression_data',sub_a,'*_mirr_mRNA.nii'))
    GABAa_mrna = ni.load(GABAa_path[0]).get_fdata()
    #finding max value for whole brain to normalise data once region has been extracted
    GABAa_max = GABAa_mrna.max()
    #Indexing each expression by the mask- shows only the receptor expression within the MRS ROI
    region_GABAa = GABAa_mrna[region_mask==1]
    #normalising region mrna values by whole brain receptors
    region_GABAa_norm = region_GABAa/GABAa_max
    #array for all normalised values for all GABAa subunits
    GABAa_region_data[:,a] = region_GABAa_norm

#convert numpy array to pandas dataframe    
GABAa_region = pd.DataFrame(data=GABAa_region_data) 
#label columns with subunit names   
GABAa_region.columns = [GABAa]
#remove any rows with 0 values
GABAa_region_removed = GABAa_region[(GABAa_region != 0).all(1)]


#GABAa Violin plots
sns.set_theme(style="whitegrid")
GABAa_plot= sns.violinplot(data=GABAa_region_removed, inner = "point")

