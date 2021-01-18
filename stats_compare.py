# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:57:08 2021

@author: lizmc
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import seaborn as sns
from scipy.stats import rankdata
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi 
from matplotlib.ticker import AutoMinorLocator

def stats_compare(data, receptors, subunits, masks):
    
    alldata_outputs = {}
    alldata_reorder = {}
    
    
    #creating blank figure for comparison plots
    out_pdf = pdf.PdfPages("regional_comparison.pdf");
    fig,axs = plt.subplots(nrows=3, ncols=2,sharex=False,sharey=False, figsize=(35, 44))
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.08, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=30)
    a=0
    b=0
    
    #loop through receptors
    for i, r in enumerate(receptors):
    #subunits list for the current receptor
        r_sub = subunits.loc[subunits[1] == r]
        aov_table = {}
        pvals = []
        fvals = []
    
        compare_data = pd.DataFrame(columns=['value','mask','subunit'])
       
        #loop subunits
        for x, sub in enumerate(r_sub[0]):
            
         
            #loop masks
            for z, m in enumerate(masks):
                #restructuring data to compare subunits between masks
                data_col = pd.DataFrame(data[m][r][:,x])
                data_col.columns = ['value']
                data_col['mask'] = m
                #subunit data for one mask
                data_col['subunit'] = sub
                #subunit data for all masks
                compare_data = pd.concat([compare_data,data_col])
                

            #one way anovas compare each subunits for all masks
            model = ols('value ~ C(mask)', data=compare_data).fit()
            aov_table[sub] = sm.stats.anova_lm(model, typ=2)
            pvals.append(aov_table[sub]['PR(>F)']['C(mask)'])
            fvals.append(aov_table[sub]['F']['C(mask)'])
        
        #corrections for multiple comparisons
        pvals_numpy = np.array(pvals).T
        corrected_pvals = multi.multipletests(pvals_numpy, alpha=0.05, method='bonferroni')
        orig_pvals = pd.DataFrame(pvals)
        orig_pvals = orig_pvals.T
        new_pvals = pd.DataFrame(corrected_pvals[1])
        new_pvals = new_pvals.T
        fvals = pd.DataFrame(fvals)
        fvals = fvals.T
        subunit_outputs = pd.concat([fvals, orig_pvals, new_pvals]) 
        subunit_outputs.columns = [r_sub[0]]   
        subunit_outputs.index = ['F','p', 'Adjusted p']

        alldata_reorder[r]=compare_data
        alldata_outputs[r]=subunit_outputs
        
        #remove 0's and make violin plots
        df_removed = compare_data[(compare_data != 0).all(1)]
        plot = sns.violinplot(x="subunit", y="value", hue="mask", data=df_removed, ax=axs[b,a])
        plot.set_title(r, fontsize=26)
        plot.tick_params(axis='both', which='both', labelsize=18)
        plot.set_xticklabels([])
        plot.set_xlabel("")
        plot.set_ylabel("")
        plot.legend(loc=1, prop={'size':20})
        minor_locator = AutoMinorLocator(2)
        plot.xaxis.set_minor_locator(minor_locator)

        plot.grid(which='minor', linestyle='-', linewidth='0.5')
    
        #adding table below figures
        cell_text = []
        tmp = alldata_outputs[r].to_numpy()
        
        for row in range(0, len(tmp)):
            cell_text.append(['%1.4f' % (value) for value in tmp[row,:]])
        rows = ['F','p', 'Adjusted p']
        columns = ['%s' % (unit) for unit in r_sub[0]]
        the_table = plot.table(cellText=cell_text, 
                  rowLabels=rows, 
                  colLabels=columns,
                  loc='bottom')
        the_table.set_fontsize(16)
        the_table.scale(1,2)
        
      
        
        a=a+1 
        if a > 1: 
            b=b+1
            a=0
             
    out_pdf.savefig()   
    out_pdf.close()
       
       
    
    return  subunit_outputs