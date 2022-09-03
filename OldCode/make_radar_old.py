
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 14:47:36 2021

@author: lizmc
"""
import matplotlib.pyplot as plt
import pandas as pd
from math import pi


def make_radar(data, receptor, mask_names):
    df_t = pd.DataFrame(data) 
    df = df_t.T
    N = len(receptor)
      
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
     
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    plt.xticks(angles[:-1], receptor)
    ax.set_rlabel_position(0)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0,","0.2", "0.4", "0.6", "0.8", "1"], color="grey", size=7)
    plt.ylim(0,1)
    
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