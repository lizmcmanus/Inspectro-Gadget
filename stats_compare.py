# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:57:08 2021

@author: lizmc
"""
import numpy 
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

from scipy.stats import rankdata


def stats_compare(data, receptors, subunits, masks):

    alldata_outputs = {}
    #loop through receptors
    for i, r in enumerate(receptors):
    #subunits list for the current receptor
        r_sub = subunits.loc[subunits[1] == r]
        aov_table = {}
        pvals = []
        fvals = []
    
        #loop subunits
        for x, sub in enumerate(r_sub[0]):
            compare_data = pd.DataFrame(columns=['value','mask','subunit'])
            
            #loop masks
            for z, m in enumerate(masks):
                
                m = m.replace(".nii.gz","")
                data_col = pd.DataFrame(data[m][r][:,0])
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
        pvals = pd.DataFrame(pvals)
        pvals = pvals.T
        fdr_correct = fdr(pvals)
        fvals = pd.DataFrame(fvals)
        fvals = fvals.T
        subunit_outputs = pd.concat([fvals, pvals, fdr_correct]) 
        subunit_outputs.columns = [r_sub[0]]   
        subunit_outputs.index = ['F','p', 'FDR p']

        
        alldata_outputs[r]=subunit_outputs
                
    return alldata_outputs


          
#fdr correction function       
def fdr(p_vals):
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    
    return fdr