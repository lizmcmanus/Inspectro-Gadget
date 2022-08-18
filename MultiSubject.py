# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 21:03:16 2021

@author: lizmc

GoGoGadget Multi Participants
"""

import os
import sys
import numpy as np
import pandas as pd
import nibabel as ni
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
from extract_mrna import extract_mrna
from ex_in_multi import ex_in_multi
from make_radar import make_radar
from make_violin import make_violin
import textwrap as twp

# Project directory
#data_dir = 'E://Taiwan/Inspectro-Gadget/'
data_dir = os.getcwd()

# Directory for outputs
output_dir = os.path.join(data_dir,'expression')

#variable for all participant masks and receptors
GABA_Glu = pd.read_csv(os.path.join(data_dir,'GABA_Glu.tsv'), delimiter='\t', header=None)
receptors = pd.read_csv(os.path.join(data_dir,'receptors.tsv'), delimiter='\t', header=None)
subjects = pd.read_csv(os.path.join(data_dir,'subjects.tsv'), delimiter='\t', header=None)
#subject = np.array(subjects)
sub_id = [sub.replace(".nii.gz","") for i,sub in enumerate(subjects[0])]
sub_id = [sub2.replace("ACC_in_MNI_","") for sub2 in sub_id]
sub_id = pd.DataFrame(sub_id)
sub_id_mean = sub_id.append(['Mean'], ignore_index=True)



#creating variable for arrays and pdfs
sub_GABAGlu_data = {}
all_GABAGlu_data = {}
sub_receptors_data = np.zeros([len(receptors),len(subjects)])
radar_pdf = pdf.PdfPages("MultiSubject_radar_plots.pdf");
violin_pdf = pdf.PdfPages("MultiSubject_GABA-GLu_Violin_plots.pdf");
ExIn_pdf = pdf.PdfPages("Excitationinhibition_Vales.pdf");
ExIn = {}
ExIn_all = {}



####extracting data for all subjects
#loop through participants list to extract data
for m, mask in enumerate(subjects[0]):
    mask_filename = mask
    # Load MRS region
    region_mask = ni.load(os.path.join(data_dir, mask_filename)).get_fdata()
    mask_size = region_mask[region_mask==1]
    mask_size = len(mask_size)
    mask_filename = mask_filename.replace(".nii.gz","")


### Making array of data for each subject for all GABA/Glu subunits  
    subject_subtype_data = np.zeros([mask_size,len(GABA_Glu)])
    recep_removed = pd.DataFrame([])
   
    for f, subtype in enumerate(GABA_Glu[0]):
        col_df=pd.DataFrame([])
        subtype_path = os.path.join(data_dir,'mRNA_images',subtype+'_mirr_mRNA.nii')
        subject_subtype_data[:,f]= extract_mrna(subtype_path, region_mask)
     
        
     #removed all 0 values (that aren't actually 0's due to sigmoid normalisation) 
        col = subject_subtype_data[:,f]
        col_rem = col[col > 0.1]
        col_df[subtype] = col_rem
        #col_df.rename(columns=subtype)
        recep_removed = pd.concat([recep_removed, col_df], axis=1)
       # recep_array= np.array(recep_removed)
     
     
    #make a dictionary for each subject mask and the GABAGlu receptors arrays 
    sub_GABAGlu_data[mask] = recep_removed
    
    
    
    
    subject_receptor_data = np.zeros([mask_size,len(receptors)])   
    allrecep_removed = pd.DataFrame([])
####RADAR PLOTS FOR OTHER RECEPTORS    
    for r, receptor in enumerate(receptors[0]):
         receptors_path = os.path.join(data_dir,'mRNA_images',receptor+'_mirr_mRNA.nii')
         subject_receptor_data[:,r] = extract_mrna(receptors_path, region_mask)
         
         #removed all 0 values (that aren't actually 0's due to sigmoid normalisation) 
         col_r = subject_receptor_data[:,r]
         col_removed = col_r[col_r > 0.1]
         recep_col_df = pd.DataFrame(col_removed)
         allrecep_removed = pd.concat([allrecep_removed, recep_col_df], axis=1)
        
         #average all the subunits in each group to give 1 value per group for each voxel in the mask
         r_mean_data = allrecep_removed.mean(axis=0)
         
     #make an array for means of subunits in each subject mask  
    sub_receptors_data[:,m] = r_mean_data   

sub_ids = sub_id.values.tolist()
make_radar(sub_receptors_data, receptors[0], sub_ids)
radar_pdf.savefig()
radar_pdf.close()
    
##### VIOLIN PLOTS #####

  #creating blank figure for subunit subplots
  #out_pdf = pdf.PdfPages(mask_filename+"_receptor_expressions.pdf");
fig,axs = plt.subplots(nrows=2, ncols=3, sharex=False,sharey=False, figsize=(10, 7), linewidth=0.01)
plt.tick_params(bottom=False, top=False, left=False, right=False)
fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=8)
fig.tight_layout(pad = 4.0)
x = 0
y = 0


#this loop creates arrays for each receptor subunit with all participants int it so this cna be used to make the violins    
for re, GABAGlu_recep in enumerate(GABA_Glu[0]):
    recep_sub_data = []
    recep_allsub = pd.DataFrame([])
    
    for n, n_mask in enumerate(subjects[0]):
        sub_data= np.array(sub_GABAGlu_data[n_mask])
        recep_sub_data = sub_data[:,re]
        recep_sub_data = pd.DataFrame(recep_sub_data)
        recep_allsub = pd.concat([recep_allsub, recep_sub_data], axis=1)
        
    mean = recep_allsub.mean(axis=1)
    recep_withmean = pd.concat([recep_allsub, mean], axis=1)
   
    plot = make_violin(recep_withmean, sub_id_mean, GABAGlu_recep, axs[y,x])
    plot.set_xticklabels('')
    plot.set_xticks([])
    #adding values table below subplots
    
    sub_mean = recep_withmean.mean(axis=0)
    
    #calcuating distance from mean for each subject
    dfm = pd.DataFrame(sub_mean - sub_mean['Mean',])
    dfm_t = dfm.T
    dfm_t = np.array (dfm_t) ####### THIS IS STUPIDLY INELEGANT, BUT IT'S THE ONLY WAY I COULD GET THE DAMN ARRAY TO TRANSPOSE
    dfm_t = np.round(dfm_t, 2)

    row_lab = [twp.fill('distance from mean', 10)]
    columns = sub_id_mean[0]
    sub_table = plot.table(cellText=dfm_t,
          rowLabels= row_lab,
          ####height = 2,
          colLabels=columns,
          loc='bottom')
    sub_table.set_fontsize(6)  


    x += 1
    if x > 2:
        y += 1
        x = 0
        
    if y > 1:
        y = 0
        violin_pdf.savefig()
        for x1 in range(0,3):
          for y1 in range(0,2):
              axs[y1,x1].clear()
              
    
violin_pdf.savefig()
violin_pdf.close()

#calculating excitation/inhibition ratios for each subject mask
ExIn_all = np.zeros([1,len(subjects[0])])
for ei, mask_ei in enumerate(subjects[0]):
    data = sub_GABAGlu_data[mask_ei]
    ExIn = ex_in_multi(data)
    ExIn_all[:,ei] = ExIn

#plotting e/i ratio's in table and saving as pdf
ExIn_all = np.round(ExIn_all, 3)
title = [twp.fill('Excitation/Inhibition Ratio', 22)]
fig, ax = plt.subplots() 
ax.set_axis_off()
cols = sub_id[0] 
table = ax.table( 
    cellText = ExIn_all,  
    rowLabels = title,
    colLabels = cols,  
    loc ='center')     
fig.tight_layout(pad = 4.0)#still having trouble with massive row label!    
ax.set_title('Multi-Subjet Excitation/Inhibition Ratio''s', fontweight ="bold", loc="centre")
   
ExIn_pdf.savefig()
ExIn_pdf.close()
         
    




