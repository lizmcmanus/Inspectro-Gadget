# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Tool to show the receptor expressions within an MRS region

MRS regions masks must be (in MNI space?))

Arguments need to be given for region mask names when using command line
GABA and Glutamate regions must ine in a .tsv file "receptors.tsv"
"""
import os
import numpy as np
import pandas as pd
import nibabel as ni
from extract_mrna import extract_mrna
from make_violin import make_violin
from excitation_inhibition import ex_in
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
import sys

# Project directory
data_dir = 'E://Taiwan/Inspectro-Gadget/'

# Directory for outputs
output_dir = os.path.join(data_dir,'expression')
#variable for all masks and receptors
all_region_data = {}
# Loop through arguments for MRS regions
for argi in range (1, len(sys.argv)):
    mask_filename = sys.argv[argi]
    # Load MRS region
    region_mask = ni.load(os.path.join(data_dir, mask_filename)).get_fdata()
    mask_size = region_mask[region_mask==1]
    mask_size = len(mask_size)
    mask_filename = mask_filename.replace(".nii.gz","")
    # Load in the list of receptors and subunits
    receptors_list = pd.read_csv(data_dir+'receptors.tsv', delimiter='\t', header=None)
    #extract unique receptors
    receptors = receptors_list[1].unique()
    
    #creating blank figure for subunit subplots
    out_pdf = pdf.PdfPages(mask_filename+"_receptor_expressions.pdf");
    fig,axs = plt.subplots(nrows=3, ncols=2,sharex=False,sharey=False, figsize=(26, 30))
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.08, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=25)
    x=0
    y=0
    
    #variable to hold all subunits data
    all_subunit_data = {} 
    #variable for excitation inhibition scores per region
    ExIn = {}
    
    # loop thorugh each receptor type
    for i, receptor_t in enumerate(receptors):
        
        #makes variables for subunits of each specific receptor type
        receptor_type = receptors_list.loc[receptors_list[1] == receptor_t]
        receptor_region_data = np.zeros([mask_size,len(receptor_type[0])])
        
        #loop through subunits within each receptor type
        for a, sub_a in enumerate(receptor_type[0]):
            #call function to extract mrna data within each region
            receptor_region_data[:,a] = extract_mrna(sub_a, region_mask)
       
        #all subunits for each receptor type stored in this variable
        all_subunit_data[receptor_t]=receptor_region_data
            
        #function to make violin plots per receptor and saves all as one pdf
        make_violin(receptor_region_data, receptor_type, receptor_t, axs[y,x])
        x=x+1 
        if x > 1: 
            y=y+1
            x=0
            
    #calculating excitation inhibition scores for each MRS region
    ExIn = ex_in(all_subunit_data)
    #adding ratio to the bottom of the figure
    fig.text(0.5, 0.07, 'Regional Excitation/Inhibition Ratio = {0}'.format(ExIn), va='center', ha='center', fontsize=32)
    
    #saving pdf and dataframes for each mask files receptor data        
    all_region_data[mask_filename]=all_subunit_data 
    out_pdf.savefig()   
    out_pdf.close()
    
  
   
 
 