# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 20:32:00 2020

@author: lizmc
"""
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def make_violin(data, name, receptor):
    plt.clf();
    n=len(name[0])
    subunits=name[0]
    x=list(range(0,n))
    df = pd.DataFrame(data) 
    #label columns with subunit names   
    df.columns = [name[0]]
    subunits = name[0]
    #remove any rows with 0 values (may be that voxels were white matter so have no receptors)
    df_missing_removed = df[(df != 0).all(1)]
    
    # Violin plots of each subunit per receptor
    sns.set_theme(style="whitegrid")
    receptor_plot = sns.violinplot(data=df_missing_removed, inner = "box")
    plt.title(receptor)
    plt.ylabel("Normalised mRNA Expression Value")
    plt.xticks(x,subunits)
    plt.xticks(fontsize=7, rotation = 45) 
    
