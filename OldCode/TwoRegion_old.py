# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 20:48:04 2021

@author: lizmc

GoGoGadget multi region analysis
"""

import os
import sys
import numpy as np
import pandas as pd
import nibabel as ni
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
from extract_mrna import extract_mrna
from make_violin import make_violin
from excitation_inhibition import ex_in
from stats_compare import stats_compare
from make_radar import make_radar



# Project directory
#data_dir = 'E://Taiwan/Inspectro-Gadget/'
data_dir = os.getcwd()

# Directory for outputs
output_dir = os.path.join(data_dir,'expression')
#variable for all masks and receptors
all_region_data = {}
receptors_list = pd.read_csv(os.path.join(data_dir,'GroupedReceptors.tsv'), delimiter='\t', header=None)
subtype = pd.read_csv(os.path.join(data_dir,'subtype_receptors.tsv'), delimiter='\t', header=None)
#extract unique receptors
receptors = receptors_list[1].unique()
masks_receptor_data = np.zeros([len(receptors),2])
masks_subtype_data = np.zeros([len(subtype),2])
radar_pdf = pdf.PdfPages("radar_plots.pdf");

masks = sys.argv[1:]
# Loop through arguments for MRS regions
for m, mask in enumerate(masks):
    mask_filename = mask
    # Load MRS region
    region_mask = ni.load(os.path.join(data_dir, mask_filename)).get_fdata()
    mask_size = region_mask[region_mask==1]
    mask_size = len(mask_size)
    mask_filename = mask_filename.replace(".nii.gz","")

    
    #creating blank figure for subunit subplots
    out_pdf = pdf.PdfPages(mask_filename+"_receptor_expressions.pdf");
    fig,axs = plt.subplots(nrows=2, ncols=3, sharex=False,sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    x = 0
    y = 0

    #variable to hold all subunits data
    all_subunit_data = {}
    all_subtype_data = {}
    #variable for excitation inhibition scores per region
    
    ExIn = {}

    # loop thorugh each receptor type for violin plots
    for i, receptor_t in enumerate(receptors):

        #makes variables for subunits of each specific receptor type
        receptor_type = receptors_list.loc[receptors_list[1] == receptor_t]
        receptor_region_data = np.zeros([mask_size,len(receptor_type[0])])

        #loop through subunits within each receptor type
        for a, sub_a in enumerate(receptor_type[0]):
            #call function to extract mrna data within each region
            subunit_path = os.path.join(data_dir,'mRNA_images',sub_a+'_mirr_mRNA.nii')
            receptor_region_data[:,a] = extract_mrna(subunit_path, region_mask)

        #all subunits for each receptor type stored in this variable
        all_subunit_data[receptor_t]=receptor_region_data

        #function to make violin plots per receptor and saves all as one pdf
        make_violin(receptor_region_data, receptor_type, receptor_t, axs[y,x])
        x += 1
        if x > 2:
            y += 1
            x = 0
            
        if y > 1:
            y = 0
            out_pdf.savefig()
            for x1 in range(0,3):
              for y1 in range(0,2):
                  axs[y1,x1].clear()
        

    
    #saving pdf and dataframes for each mask files receptor data
    all_region_data[mask_filename]=all_subunit_data

    #calculating excitation inhibition scores for each MRS region
    ExIn = ex_in(all_subunit_data)
    #adding ratio to the bottom of the figure
    fig.text(0.5, 0.04, 'Regional Excitation/Inhibition Ratio = {:.2f}'.format(ExIn), va='center', ha='center', fontsize=12)

    out_pdf.savefig()
    out_pdf.close()
 
 ###radar plot stuff   
    all_subtype_data = np.zeros([mask_size,len(subtype)])
 
    subtype_region_data = np.zeros([mask_size,len(subtype)])
    
    for f, sub in enumerate(subtype[0]):
        subtype_path = os.path.join(data_dir,'mRNA_images',sub+'_mirr_mRNA.nii')
        subtype_region_data[:,f] = extract_mrna(subtype_path, region_mask)
        #average all the subunits in each group to give 1 value per group for each voxel in the mask
        mean_data = subtype_region_data.mean(axis=0)
    
        #create data frame with mean data for all grouped receptors
    all_subtype_data=mean_data
   
    #make an array for means of receptors in each mask  
    masks_subtype_data[:,m] = all_subtype_data
 
mask_names = [m.replace('.nii.gz', '') for m in masks]

make_radar(masks_subtype_data, subtype[0], mask_names)
radar_pdf.savefig()
radar_pdf.close()

#stats to compare subunit distributions between regions
stats_out = stats_compare(all_region_data, receptors, receptors_list, mask_names)