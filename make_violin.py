# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 20:32:00 2020

@author: lizmc
"""
import seaborn as sns
import pandas as pd
import numpy as np


def make_violin(data, name, receptor, fig):
    n=len(name[0])
    subunits=name[0]
    df = pd.DataFrame(data) 
    #label columns with subunit names   
    df.columns = [name[0]]
    subunits = name[0]
    #remove any rows with 0 values (may be that voxels were white matter so have no receptors)
   
    #due to new mrna extraction normalisation, 0's aren't 0, they're just very small numbers
    #df_missing_removed = df[(df >0.001).all(1)] # this is removing any row that has a missing value, not just the missing values in a column!!!
   
    
    # Violin plots of each subunit per receptor
   # sns.set_theme(style="whitegrid")
    plot = sns.violinplot(data=df, inner = "box", ax=fig, linewidth=0.1, grid_linewidth=1)
    fig.set_title(receptor, fontsize=6)
    fig.set_xticks(np.arange(0, n))
    fig.set_xticklabels(subunits, fontsize=6, rotation = 45)
    fig.tick_params(axis='y', which='major', labelsize=6)
    
    return(plot)
  

