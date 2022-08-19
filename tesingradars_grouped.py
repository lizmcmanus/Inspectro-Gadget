# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 15:38:05 2021

@author: lizmc
"""
import os
import sys
import numpy as np
import pandas as pd
import nibabel as ni
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
from extract_mrna import extract_mrna
from make_radar import make_radar
from math import pi

data_dir = os.getcwd()

masks = sys.argv[1:]
gr_receptors_list = pd.read_csv(os.path.join(data_dir,'GroupedReceptors.tsv'), delimiter='\t', header=None)
receptors = gr_receptors_list[1].unique()
    
masks_receptor_data = np.zeros([len(receptors),2])

# Loop through arguments for MRS regions
for m, mask in enumerate(masks):
    mask_filename = mask
    # Load MRS region
    region_mask = ni.load(os.path.join(data_dir, mask_filename)).get_fdata()
    mask_size = region_mask[region_mask==1]
    mask_size = len(mask_size)
    mask_filename = mask_filename.replace(".nii.gz","")
    # Load in the list of receptors and subunits

    radar_pdf = pdf.PdfPages("radar_plots.pdf");
    
    all_grouped_receptor_data = np.zeros([mask_size,len(receptors)])
    for gr, gr_t in enumerate(receptors):
        
         #makes variables for subunits of each specific receptor type
        grouped_receptors = gr_receptors_list.loc[gr_receptors_list[1] == gr_t]
        grouped_receptor_region_data = np.zeros([mask_size,len(grouped_receptors [0])])
        
        for f, sub_grouped in enumerate(grouped_receptors[0]):
            gr_subunit_path = os.path.join(data_dir, 'inspectro_gadget/data/mRNA_images', sub_grouped + '_mirr_mRNA.nii')
            grouped_receptor_region_data[:,f] = extract_mrna(gr_subunit_path, region_mask)
            #average all the subunits in each group to give 1 value per group for each voxel in the mask
            group_mean_data = grouped_receptor_region_data.mean(axis=1)
        
            #create data frame with mean data for all grouped receptors
        all_grouped_receptor_data[:,gr]=group_mean_data
        
    all_receptor_data_removed = all_grouped_receptor_data[(all_grouped_receptor_data != 0).all(1)]
    #mean of all recptors in the mask
    gr_data_means = all_receptor_data_removed.mean(axis=0)
   
    #make an array for means of receptors in each mask  
    masks_receptor_data[:,m] = gr_data_means 

mask_names = [m.replace('.nii.gz', '') for m in masks]
    
df_t = pd.DataFrame(masks_receptor_data) 
df = df_t.T
N = len(receptors)
  
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]
 
fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)
plt.xticks(angles[:-1], receptors)
ax.set_rlabel_position(0)
plt.yticks([4.5, 4.75, 5, 5.25, 5.5], ["4.5","4.75","5","5.25","5"], color="grey", size=7)
plt.ylim(4.5,5.5)

# mask 1
values=df.loc[0].values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label=mask_names[0])
ax.fill(angles, values, 'b', alpha=0.1)
 
# mask 2
values=df.loc[1].values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label=mask_names[1])
ax.fill(angles, values, 'r', alpha=0.1)
 
# Add legend
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))


radar_pdf.savefig()
radar_pdf.close()
   

       



       

   
   
