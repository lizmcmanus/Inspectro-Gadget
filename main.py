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

# Project directory
data_dir = 'E://Taiwan/Inspectro-Gadget/'

# Directory for outputs
output_dir = os.path.join(data_dir,'expression')

# Create outputs directory if missing
#if not os.path.isdir(output_dir):
#  os.mkdir(output_dir)

# Load MRS region
region_mask = ni.load(os.path.join(data_dir,'MRS_region.nii.gz')).get_fdata()

# Load in the list of receptors
receptors = pd.read_csv(data_dir+'receptors.tsv', delimiter='\t')
n_receptors = len(receptors)

#creating dataframe for all receptor data
all_receptor_data = np.zeros(np.hstack((n_receptors,region_mask.shape)))
print(all_receptor_data)

# Load data for each recetpors
for i,sub in enumerate(receptors):
    receptor_path = glob.glob(os.path.join(data_dir,'mRNA_Expression_data',sub,'*_mirr_mRNA.nii'))
    mrna = ni.load(receptor_path[0]).get_fdata()
    #Indexing each expression by the mask
    #shows only the receptor expression within the MRS ROI
    region_values = mrna[region_mask==1]
    all_receptor_data[i,:,:,:] = mrna
