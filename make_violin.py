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
    df_missing_removed = df[(df != 0).all(1)]
    
    # Violin plots of each subunit per receptor
   # sns.set_theme(style="whitegrid")
    sns.violinplot(data=df_missing_removed, inner = "box", ax=fig)
    fig.set_title(receptor, fontsize=22)
    fig.set_xticks(np.arange(0, n))
    fig.set_xticklabels(subunits, fontsize=14, rotation = 45)
    fig.tick_params(axis='y', which='major', labelsize=14)

