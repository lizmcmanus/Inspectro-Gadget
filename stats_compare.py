# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:57:08 2021

@author: lizmc
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi
from matplotlib.ticker import AutoMinorLocator
from scipy.stats import trim_mean, norm, mstats
from joblib import Parallel, delayed
from scipy import stats
from sklearn import preprocessing
import textwrap as twp


def stats_compare(data, receptors, subunits, mask_names,n_samples=10000):
    print('Running region comparisons for:')

    alldata_outputs = {}
    alldata_reorder = {}


    #creating blank figure for comparison plots
    out_pdf = pdf.PdfPages("regional_comparison.pdf");
    fig,axs = plt.subplots(nrows=2, ncols=3, sharex=False,sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015,0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    a=0
    b=0

    #loop through receptors
    for i, r in enumerate(receptors):

    #subunits list for the current receptor
        r_sub = subunits.loc[subunits[1] == r]
        aov_table = {}

        dvals = []
        dcis_up = []
        dcis_low = []
        pctDifs = []
        KStest_all = []

        compare_data = pd.DataFrame(columns=['value','mask','subunit'])

        alpha = 0.05/len(r_sub)
        #loop subunits
        for x, sub in enumerate(r_sub[0]):
           # print(f'{r} - {sub}')
            # Get data for each mask
            mask1_data = data[mask_names[0]][r][:,x]
            mask2_data = data[mask_names[1]][r][:,x]
            
            # Remove zeros
            # Do we want to do this though? The zeros skew the distribution so
            # outliers of non-zero voxels may not be picked up but at the same
            # time the zeros are data points from the voxel
            mask1_data = mask1_data[mask1_data != 0]
            mask2_data = mask2_data[mask2_data != 0]
            
            #mean centering data for comparison of distributions
            mc_mask1 = mask1_data - np.mean(mask1_data)
            mc_mask2 = mask2_data - np.mean(mask2_data)
            #comparing distribution
            KStest = stats.kstest(mc_mask1, mc_mask2,  alternative ='two-sided', mode = 'auto')
            KStest = KStest[0]
            
            
            # Compare subunit values
            pctDif,D,Dci_low,Dci_up = bootstrap_diff(mask1_data,mask2_data,n_samples=n_samples,alpha=alpha)

            dvals.append(D)
            dcis_up.append(Dci_up)
            dcis_low.append(Dci_low)
            pctDifs.append(pctDif)
            KStest_all.append(KStest)
            
           
            ##### added this back in so the data is used for the violins 
            for z, m in enumerate(mask_names):
                #restructuring data to compare subunits between masks
                data_col = pd.DataFrame(data[m][r][:,x])
                data_col.columns = ['value']
                data_col['mask'] = m
                #subunit data for one mask
                data_col['subunit'] = sub
                #subunit data for all masks
                compare_data = pd.concat([compare_data,data_col])
         
            
         #### string cohens D and CI's 
            dval_row=[]
            format="{d:.2f} ({dci_low: .2f} - {dci_up: .2f})"

            for d in range(0,len(dvals)):
                dval_row.append(format.format(d=dvals[d], dci_low=dcis_low[d], dci_up=dcis_up[d]))
            dval_row=pd.DataFrame(dval_row)
          
               
    
    
    ######################################

        #corrections for multiple comparisons
        # pvals_numpy = np.array(pvals).T
        # corrected_pvals = multi.multipletests(pvals_numpy, alpha=0.05, method='bonferroni')
        # orig_pvals = pd.DataFrame(pvals).T
        # new_pvals = pd.DataFrame(corrected_pvals[1]).T
        #dvals = pd.DataFrame(dvals)
        #pctDifs_df = pd.DataFrame(pctDifs)
        pctDifs_row = []
        for pc in range(0, len(pctDifs)):
            pctDifs_row.append("%.2f" % pctDifs[pc])
        pctDifs_df = pd.DataFrame(pctDifs_row)
        #dcis_low_df = pd.DataFrame(dcis_low)
       # dcis_up_df = pd.DataFrame(dcis_up)
        KS_row = []
        for ks in range(0, len(KStest_all)):
            KS_row.append("%.2f" % KStest_all[ks])
        KStest_df = pd.DataFrame(KS_row)
        
        
      
        ##################################
        
        
        
        subunit_outputs = pd.concat([pctDifs_df, dval_row, KStest_df], axis=1)
        subunit_outputs = subunit_outputs.T
        subunit_outputs.columns = [r_sub[0]]
        
        

        alldata_reorder[r]=compare_data
        alldata_outputs[r]=subunit_outputs

        #remove 0's and make violin plots
        #df_removed = compare_data[(compare_data != 0).all(1)]
        df_removed = compare_data[(compare_data['value'] > 0.1)]

        
        sns.set_palette("hls", 8)
        #axs[b,a].clear()
        plot = sns.violinplot(x="subunit", y="value", hue="mask", data=df_removed, ax=axs[b,a], linewidth=0.1, grid_linewidth=1)
        plot.set_title(r, fontsize=6)
        plot.tick_params(axis='y', which='both', labelsize=4, width=0.5)
        plot.set_xticklabels([])
        plot.set_xticks([])
        plot.set_xlabel("")
        plot.set_ylabel("")
        plot.legend(loc=0, prop={'size':6})
        minor_locator = AutoMinorLocator(2)
        plot.xaxis.set_minor_locator(minor_locator)
        plot.grid(which='minor', linestyle='-', linewidth='0.01')
        fig.tight_layout(pad=3.0)

        #adding table below figures
        cell_text = []
        tmp = subunit_outputs.to_numpy()

        #cell_text[[f'{r}']] ###WTF does this do?!

        for row in range(0, len(tmp)):
            cell_text.append(['%s' % (value) for value in tmp[row,:]])
        rows = ['% difference','Cohen\'s d', "KS_D"]
        columns = ['%s' % (unit) for unit in r_sub[0]]
        the_table = plot.table(cellText=cell_text,
              rowLabels=rows,
              colLabels=columns,
              loc='bottom')
        the_table.set_fontsize(8)
        the_table.scale(1,0.7)
        
       ### subplot counters     
        a += 1 
        if a > 2:
            b += 1
            a = 0
            
        if b > 1:
            b = 0
            out_pdf.savefig()
            for a1 in range(0,3):
               for b1 in range(0,2):
                   axs[b1,a1].clear()
        
       
            
    out_pdf.savefig()
    out_pdf.close()


    return  subunit_outputs

def calc_cohend(g1,g2):
    """
    Calculate Cohen's d for a comparison of two dependent groups.

    g1: First group data
    g2: Second group data

    """
    m1 = np.mean(g1)
    m2 = np.mean(g2)
    sd = np.std(g1)+np.std(g2)/2
    #sd = np.std(mstats.winsorize(g1,limits=(0.1,0.1)))+np.std(mstats.winsorize(g2,limits=(0.1,0.1)))/2
    return((m1-m2)/sd)

def get_boot_d(g1,g2):
    """
    Calculate Cohen's d for a single pair of resampled data.

    g1: First group data
    g2: Second group data

    """
    idx1 = np.random.randint(0,g1.shape[0],g1.shape[0])
    idx2 = np.random.randint(0,g2.shape[0],g2.shape[0])
    return(calc_cohend(g1[idx1],g2[idx2]))

def get_jack_d(g1,g2,i):
    """
    Calculate Cohen's d for a single pair of jackknife data.

    g1: First group data
    g2: Second group data
    i: iteration number

    """
    j1 = np.delete(g1,i,axis=0)
    j2 = np.delete(g2,i,axis=0)
    return(calc_cohend(j1,j2))

def cohens_d(g1,g2,n_samples=10000,alpha=0.05):
    """
    Calculate Cohen's d plus adjusted confidence interval for a comparison
    of two dependent groups.

    g1: First group data
    g2: Second group data
    n_samples: number of permutations
    alpha: confidence interval width

    """
    alphas = np.array([alpha/2,1-alpha/2])
    orig = calc_cohend(g1,g2)
    boot = sorted(Parallel(n_jobs=-2,verbose=0)(delayed(get_boot_d)(g1,g2)for i in range(n_samples)))
    z0 = norm.ppf((1.0*np.sum(boot < orig, axis=0))/n_samples)
    jstat = Parallel(n_jobs=-2,verbose=0)(delayed(get_jack_d)(g1,g2,i)for i in range(min((len(g1),len(g2)))))
    jmean = np.mean(jstat,axis=0)
    a = np.sum( (jmean - jstat)**3, axis=0 ) / ( 6.0 * np.sum( (jmean - jstat)**2, axis=0)**1.5 )
    zs = z0 + norm.ppf(alphas).reshape(alphas.shape+(1,)*z0.ndim)
    avals = norm.cdf(z0 + zs/(1-a*zs))
    nvals = np.around((n_samples-1)*avals).astype('int')
    return(orig,[boot[nvals[0]],boot[nvals[1]]])

def mad_median(x):
    """
    Remove outliers based upon median absolute deviation. Threshold is set
    by chi squared distribution with two degrees of freedom.

    x: data to remove outliers from
    """
    mad = np.median(abs(x-np.median(x)))
    mm = (x-np.median(x))/(mad/0.6745)
    return(x[mm<2.24])

def perm_median(g1,g2):
    """
    Calculate median for single pair of shuffled data.

    g1: First group data
    g2: Second group data
    """
    len1 = len(g1)
    comb_gs = np.hstack((g1,g2))
    np.random.shuffle(comb_gs)
    return(np.median(comb_gs[:len1])-np.median(comb_gs[len1:]))

def bootstrap_diff(g1,g2,n_samples=10000,alpha=0.05):
    """
    Run robust comparison of two independent samples.

    Removes outliers based on median absolute deviation. Tests significance of
    difference in medians via permutation test. Calculates Cohen's d plus
    confidence interval of this. Returns these plus percentage difference between
    samples.


    g1: First group data
    g2: Second group data

    return: percentage difference, Cohen's d, Cohens's d confidence interval,
            p-value

    """
    D,Dci = cohens_d(g1,g2,alpha=alpha)
    g1_mad = mad_median(g1)
    g2_mad = mad_median(g2)
    mDif = np.median(g1_mad)-np.median(g2_mad)
    pctDif = mDif/np.median(g2_mad)*100
    mPerm = np.array(Parallel(n_jobs=-2,verbose=0)(delayed(perm_median)(g1_mad,g2_mad)for i in range(n_samples)))
    # p = len(np.where(np.abs(mPerm)>=np.abs(mDif))[0])/float(n_samples)
    return(pctDif,D,Dci[0],Dci[1])
