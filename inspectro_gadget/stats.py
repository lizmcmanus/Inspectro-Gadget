#!/usr/bin/env python3

"""
Resources for handling statistical operations and similar.
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.backends.backend_pdf as pdf
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi
from matplotlib.ticker import AutoMinorLocator
from scipy.stats import trim_mean, norm, mstats, kstest
from joblib import Parallel, delayed
from sklearn import preprocessing


def compare_regions(subunit_data, receptor_list, n_samples=5000):
    """
    Compare receptor subunit expression between two regions.

    Parameters
    ----------
    subunit_data: dict
        Dictionary with the receptor expression values for the two regions
    receptor_list: dataframe
        Dataframe containing the list of receptor subunits and the receptors they compose. Analysis is grouped by
        receptor type
    n_samples: int
        Number of permutations for calculating confidence intervals

    Returns
    -------


    """
    region_names = list(subunit_data.keys())
    subunit_d_vals = {}
    subunit_d_cis = {}
    subunit_ks_vals = {}
    subunit_pct_diff = {}
    # Loop through receptor types
    for receptor in np.unique(receptor_list['grouping']):
        # Modify alpha for number of comparisons
        alpha = 0.05/receptor_list[receptor_list.grouping == receptor].shape[0]
        for subunit in receptor_list[receptor_list.grouping == receptor]['subunit']:
            # Get data and remove NaNs
            region_one = subunit_data[region_names[0]][subunit].values
            region_two = subunit_data[region_names[1]][subunit].values
            region_one = region_one[~np.isnan(region_one)]
            region_two = region_two[~np.isnan(region_two)]
            # Compare subunit values
            subunit_pct_diff[subunit], subunit_d_vals[subunit], subunit_d_cis[subunit] = bootstrap_diff(region_one,
                                                                                                        region_two,
                                                                                                        n_samples=n_samples,
                                                                                                        alpha=alpha)
            # Mean centre for KS test
            region_one = region_one - np.mean(region_one)
            region_two = region_two - np.mean(region_two)
            # Apply KS test
            subunit_ks_vals[subunit] = kstest(region_one, region_two,  alternative='two-sided', mode='auto').statistic
    return subunit_d_vals, subunit_d_cis, subunit_pct_diff, subunit_ks_vals


def calc_cohend(g1, g2):
    """
    Calculate Cohen's d for a comparison of two dependent groups.

    g1: First group data
    g2: Second group data

    """
    m1 = np.mean(g1)
    m2 = np.mean(g2)
    sd = np.std(g1) + np.std(g2) / 2
    #sd = np.std(mstats.winsorize(g1,limits=(0.1,0.1)))+np.std(mstats.winsorize(g2,limits=(0.1,0.1)))/2
    return (m1 - m2) / sd


def get_boot_d(g1, g2):
    """
    Calculate Cohen's d for a single pair of resampled data.

    g1: First group data
    g2: Second group data

    """
    idx1 = np.random.randint(0, g1.shape[0], g1.shape[0])
    idx2 = np.random.randint(0, g2.shape[0], g2.shape[0])
    return calc_cohend(g1[idx1], g2[idx2])


def get_jack_d(g1, g2, ii):
    """
    Calculate Cohen's d for a single pair of jackknife data.

    g1: First group data
    g2: Second group data
    i: iteration number

    """
    j1 = np.delete(g1, ii, axis=0)
    j2 = np.delete(g2, ii, axis=0)
    return calc_cohend(j1, j2)


def cohens_d(g1, g2, n_samples=10000, alpha=0.05):
    """
    Calculate Cohen's d plus adjusted confidence interval for a comparison
    of two dependent groups.

    g1: First group data
    g2: Second group data
    n_samples: number of permutations
    alpha: confidence interval width

    """
    alphas = np.array([alpha/2, 1-alpha/2])
    orig = calc_cohend(g1, g2)
    boot = sorted(Parallel(n_jobs=-2, verbose=0)(delayed(get_boot_d)(g1, g2)for ii in range(n_samples)))
    z0 = norm.ppf((1.0*np.sum(boot < orig, axis=0))/n_samples)
    jstat = Parallel(n_jobs=-2,verbose=0)(delayed(get_jack_d)(g1, g2, ii)for ii in range(min((len(g1), len(g2)))))
    jmean = np.mean(jstat, axis=0)
    a = np.sum((jmean - jstat)**3, axis=0) / (6.0 * np.sum((jmean - jstat)**2, axis=0)**1.5)
    zs = z0 + norm.ppf(alphas).reshape(alphas.shape+(1,)*z0.ndim)
    avals = norm.cdf(z0 + zs/(1-a*zs))
    nvals = np.around((n_samples-1)*avals).astype('int')
    return orig, [boot[nvals[0]], boot[nvals[1]]]


def mad_median(x):
    """
    Remove outliers based upon median absolute deviation. Threshold is set
    by chi squared distribution with two degrees of freedom.

    x: data to remove outliers from
    """
    mad = np.median(abs(x-np.median(x)))
    mm = (x-np.median(x))/(mad/0.6745)
    return x[mm < 2.24]


def perm_median(g1, g2):
    """
    Calculate median for single pair of shuffled data.

    g1: First group data
    g2: Second group data
    """
    len1 = len(g1)
    comb_gs = np.hstack((g1, g2))
    np.random.shuffle(comb_gs)
    return np.median(comb_gs[:len1])-np.median(comb_gs[len1:])


def bootstrap_diff(g1, g2, n_samples=10000, alpha=0.05):
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
    D, Dci = cohens_d(g1, g2, alpha=alpha, n_samples=n_samples)
    g1_mad = mad_median(g1)
    g2_mad = mad_median(g2)
    mDif = np.median(g1_mad)-np.median(g2_mad)
    pctDif = mDif/np.median(g2_mad)*100
    #mPerm = np.array(Parallel(n_jobs=-2,verbose=0)(delayed(perm_median)(g1_mad,g2_mad)for i in range(n_samples)))
    # p = len(np.where(np.abs(mPerm)>=np.abs(mDif))[0])/float(n_samples)
    return pctDif, D, Dci
