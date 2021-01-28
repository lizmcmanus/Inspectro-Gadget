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

def stats_compare(data, receptors, subunits, mask_names,n_samples=10000):
    print('Running region comparisons for:')

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
        print(f'{r}')

    #subunits list for the current receptor
        r_sub = subunits.loc[subunits[1] == r]
        aov_table = {}
        pvals = []
        dvals = []
        dcis = []
        pctDifs = []

        compare_data = pd.DataFrame(columns=['value','mask','subunit'])

        #loop subunits
        for x, sub in enumerate(r_sub[0]):
            # Get data for each mask
            mask1_data = data[mask_names[0]][r][:,x]
            mask2_data = data[mask_names[1]][r][:,x]

            # Compare subunit values
            pctDif,D,Dci,p = bootstrap_diff(mask1_data,mask2_data,n_samples=n_samples)

            pvals.append(p)
            dvals.append(D)
            dcis.append(Dci)
            pctDifs.append(pctDif)


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

        sns.set_palette("hls", 8)
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

        a =+ 1
        if a > 1:
            b =+ 1
            a = 0

    out_pdf.savefig()
    out_pdf.close()

    return  subunit_outputs

def calc_cohend(g1,g2):
    """
    Calculate Cohen's d for a comparison of two idependent groups.

    g1: First group data
    g2: Second group data

    """
    m1 = trim_mean(g1,0.1)
    m2 = trim_mean(g2,0.1)
    sd = np.sqrt((((len(g1)-1)*np.std(mstats.winsorize(g1,limits=(0.1,0.1)))**2)+((len(g2)-1)*np.std(mstats.winsorize(g2,limits=(0.1,0.1)))**2))/(len(g1+len(g2)-2)))
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
    of two independent groups.

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

def bootstrap_diff(g1,g2,n_samples=10000):
    """
    Run robust comparison of two independent samples.

    Removes outliers based on median absolute deviation. Tests significance of
    difference in medians via permutation test. Calculates robust Cohen's d
    (20% trimmed mean & 20% Winsorized variance) plus confidence interval of this.
    Returns these plus percentage difference between samples.


    g1: First group data
    g2: Second group data

    return: percentage difference, Cohen's d, Cohens's d confidence interval,
            p-value

    """
    D,Dci = cohens_d(g1,g2)
    g1_mad = mad_median(g1)
    g2_mad = mad_median(g2)
    mDif = np.median(g1_mad)-np.median(g2_mad)
    pctDif = mDif/np.median(g2_mad)*100
    mPerm = np.array(Parallel(n_jobs=-2,verbose=0)(delayed(perm_median)(g1_mad,g2_mad)for i in range(n_samples)))
    p = len(np.where(np.abs(mPerm)>=np.abs(mDif))[0])/n_samples
    return(pctDif,D,Dci,p)
