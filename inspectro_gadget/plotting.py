#!/usr/bin/env python3

"""
Resources for plotting MRS voxel locations and the various mRNA value plots.

"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sb
from matplotlib.ticker import AutoMinorLocator
from scipy.ndimage import center_of_mass


def plot_masks(mask_imgs, labels, bground, pdf):
    """
    Show the location of the MRS voxels on a background image.

    Parameters
    ----------
    mask_imgs: list
        List containing the MRS voxel mask images
    labels: list
        List contatining the region labels
    bground: array
        Numpy array with the background image masks should be superimposed upon
    pdf: PDFPages object
        PDF object to add figures to

    Returns
    -------
    Figure object

    """
    no_ax = len(labels)
    fig = plt.figure(tight_layout=True)
    gs = gridspec.GridSpec(no_ax, 3)
    for aa in range(no_ax):
        cmass = np.round(center_of_mass(mask_imgs[aa]), 0).astype(int)
        # Sagittal
        mask_roi = np.ma.masked_where(mask_imgs[aa][cmass[0], :, :] == 0, mask_imgs[aa][cmass[0], :, :])
        ax = fig.add_subplot(gs[aa, 0])
        ax.imshow(bground[cmass[0], :, :].T, cmap='Greys_r', origin='lower', interpolation='none')
        ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
        ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
        ax.tick_params(left=False, bottom=False)
        # Coronal
        mask_roi = np.ma.masked_where(mask_imgs[aa][:, cmass[1], :] == 0, mask_imgs[aa][:, cmass[1], :])
        ax = fig.add_subplot(gs[aa, 1])
        ax.imshow(bground[:, cmass[1], :].T, cmap='Greys_r', origin='lower', interpolation='none')
        ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
        ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
        ax.tick_params(left=False, bottom=False)
        ax.set_title(labels[aa])
        # Axial
        mask_roi = np.ma.masked_where(mask_imgs[aa][:, :, cmass[2]] == 0, mask_imgs[aa][:, :, cmass[2]])
        ax = fig.add_subplot(gs[aa, 2])
        ax.imshow(bground[:, :, cmass[2]].T, cmap='Greys_r', origin='lower', interpolation='none')
        ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
        ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
        ax.tick_params(left=False, bottom=False)
    pdf.savefig(fig)
    return pdf


def plot_overlap(overlap, labels, bground, pdf):
    """
    Show the overlap between subject MRS regions on a background image.

    Parameters
    ----------
    overlap: array
        Numpy array containing overlap image
    labels: list
        List containing region name
    bground: array
        Numpy array containing background image
    pdf: PDFPages object
        PDF object to add figures to

    Returns
    -------
    Figure object

    """
    fig = plt.figure(tight_layout=True)
    gs = gridspec.GridSpec(1, 3)
    cmass = np.round(center_of_mass(overlap), 0).astype(int)
    # Sagittal
    mask_roi = np.ma.masked_where(overlap[cmass[0], :, :] == 0, overlap[cmass[0], :, :])
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(bground[cmass[0], :, :].T, cmap='Greys_r', origin='lower', interpolation='none')
    ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
    ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
    ax.tick_params(left=False, bottom=False)
    # Coronal
    mask_roi = np.ma.masked_where(overlap[:, cmass[1], :] == 0, overlap[:, cmass[1], :])
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(bground[:, cmass[1], :].T, cmap='Greys_r', origin='lower', interpolation='none')
    ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
    ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
    ax.tick_params(left=False, bottom=False)
    ax.set_title(labels[0])
    # Axial
    mask_roi = np.ma.masked_where(overlap[:, :, cmass[2]] == 0, overlap[:, :, cmass[2]])
    ax = fig.add_subplot(gs[0, 2])
    ax.imshow(bground[:, :, cmass[2]].T, cmap='Greys_r', origin='lower', interpolation='none')
    ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
    ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
    ax.tick_params(left=False, bottom=False)
    pdf.savefig(fig)
    return pdf


def make_single_violin(ax, receptors, group):
    """
    Create a subplot with violin plots showing receptor expression for a group of subunits.

    Parameters
    ----------
    ax: matplotlib axis
        The axis to add the plot to
    receptors: dataframe
        Dataframe with receptor expression. Column names are subunit names. Includes only subunits for the receptor to
        be plotted.
    group: string
        The subunit grouping name for use in the title

    Returns
    -------
    ax: matplotlib axis
        The axis with the plots added
    """
    sb.violinplot(data=receptors, inner="box", ax=ax, linewidth=0.1, grid_linewidth=1)
    ax.set_title(group, fontsize=6)
    ax.tick_params(axis='x', which='major', labelsize=6, rotation=45)
    ax.tick_params(axis='y', which='major', labelsize=6)
    return ax


def single_region_violins(subunit_data, receptor_list, pdf):
    """

    Parameters
    ----------
    subunit_data: dataframe
        Dataframe with receptor expression. Column names are subunit names
    receptor_list: dataframe
        Dataframe containing the list of receptor subunits and the receptors they compose. Plots are grouped by
        receptor type
    pdf: PDFPages object
        PDF object to add figures to

    Returns
    -------

    """
    # First PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    for rr, receptor in enumerate(['GABAA_Alpha', 'GABAA_Beta', 'GABAA_Gamma']):
        axs[0, rr] = make_single_violin(axs[0, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    for rr, receptor in enumerate(['GABAB', 'NMDA', 'AMPA']):
        axs[1, rr] = make_single_violin(axs[1, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    pdf.savefig(fig)
    plt.close()
    # Second PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    for rr, receptor in enumerate(['mGlu(I)', 'mGlu(II)', 'mGlu(III)']):
        axs[0, rr] = make_single_violin(axs[0, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    for rr, receptor in enumerate(['Kainate', 'Dopamine', '5-HT1']):
        axs[1, rr] = make_single_violin(axs[1, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    plt.close()
    pdf.savefig(fig)
    # Third PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    for rr, receptor in enumerate(['5-HT2', '5-HT3+', 'NAalpha1']):
        axs[0, rr] = make_single_violin(axs[0, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    for rr, receptor in enumerate(['NAalpha2', 'NAbeta']):
        axs[1, rr] = make_single_violin(axs[1, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    fig.delaxes(axs[1][2])
    plt.close()
    pdf.savefig(fig)
    return pdf


def make_two_violins(ax, receptors, group, pcts, ds, ds_ci, kss):
    """
    Create a subplot with violin plots comparing receptor expression between two regions for a group of subunits.

    Plots have table beneath them with comparison statistics.

    Parameters
    ----------
    ax: matplotlib axis
        The axis to add the plot to
    receptors: dataframe
        Dataframe with receptor expression. Column names are subunit names. Includes only subunits for the receptor to
        be plotted.
    group: string
        The subunit grouping name for use in the title
    stats:

    Returns
    -------
    ax: matplotlib axis
        The axis with the plots added
    """
    # Table contents
    rows = ['% diff.', "Cohen's d", 'KS']
    cell_text = [[f'{pct:0.2f}' for pct in pcts],
                 [f'{d:0.2f} ({ds_ci[ii][0]:0.2f} {ds_ci[ii][1]:0.2f})' for ii, d in enumerate(ds)],
                 [f'{ks:0.2f}' for ks in kss]]
    columns = np.unique(receptors.subunit)
    # Plot data
    sb.violinplot(data=receptors, x='subunit', y='values', hue='region',
                  inner="box", ax=ax, linewidth=0.1, grid_linewidth=1)
    ax.set_title(group, fontsize=6)
    ax.tick_params(axis='y', which='both', labelsize=4, width=0.5)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(loc=0, prop={'size': 6})
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.grid(which='minor', linestyle='-', linewidth='0.01')
    # Add table
    the_table = table(ax, cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
    the_table.set_fontsize(8)
    the_table.scale(1, 0.7)
    return ax


def two_region_prep(subunit_data, receptor_list, receptor , subunit_pct_diff, subunit_d_vals, subunit_d_cis, subunit_ks_vals):
    # Arrange subunit data into dataframe for plotting
    subunit_exp = pd.DataFrame(columns=['values', 'subunit', 'region'])
    for region in list(subunit_data.keys()):
        for subunit in receptor_list.subunit[receptor_list.grouping == receptor].values:
            df = pd.DataFrame()
            df['values'] = subunit_data[region][subunit]
            df['subunit'] = subunit
            df['region'] = region
            subunit_exp = pd.concat((subunit_exp, df), axis=0)
    subunit_exp = subunit_exp.dropna()
    # Arrange stats
    pcts = [subunit_pct_diff[subunit] for subunit in receptor_list.subunit[receptor_list.grouping == receptor].values]
    ds = [subunit_d_vals[subunit] for subunit in receptor_list.subunit[receptor_list.grouping == receptor].values]
    ds_ci = [subunit_d_cis[subunit] for subunit in receptor_list.subunit[receptor_list.grouping == receptor].values]
    kss = [subunit_ks_vals[subunit] for subunit in receptor_list.subunit[receptor_list.grouping == receptor].values]
    return subunit_exp, receptor, pcts, ds, ds_ci, kss


def two_region_violins(subunit_data, receptor_list, pdf, subunit_pct_diff, subunit_d_vals, subunit_d_cis, subunit_ks_vals):
    # First PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    for rr, receptor in enumerate(['GABAA_Alpha', 'GABAA_Beta', 'GABAA_Gamma']):
        # Prepare data
        subunit_exp, receptor, pcts, ds, ds_ci, kss = two_region_prep(subunit_data, receptor_list, receptor,
                                                                      subunit_pct_diff, subunit_d_vals,
                                                                      subunit_d_cis, subunit_ks_vals)
        # Make plot
        axs[0, rr] = make_two_violins(axs[0, rr], subunit_exp, receptor, pcts, ds, ds_ci, kss)
    for rr, receptor in enumerate(['GABAB', 'NMDA', 'AMPA']):
        # Prepare data
        subunit_exp, receptor, pcts, ds, ds_ci, kss = two_region_prep(subunit_data, receptor_list, receptor,
                                                                      subunit_pct_diff, subunit_d_vals,
                                                                      subunit_d_cis, subunit_ks_vals)
        # Make plot
        axs[1, rr] = make_two_violins(axs[1, rr], subunit_exp, receptor, pcts, ds, ds_ci, kss)
    fig.tight_layout(pad=3.0)
    pdf.savefig(fig)
    plt.close()
    # Second PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    for rr, receptor in enumerate(['mGlu(I)', 'mGlu(II)', 'mGlu(III)']):
        # Prepare data
        subunit_exp, receptor, pcts, ds, ds_ci, kss = two_region_prep(subunit_data, receptor_list, receptor,
                                                                      subunit_pct_diff, subunit_d_vals,
                                                                      subunit_d_cis, subunit_ks_vals)
        # Make plot
        axs[0, rr] = make_two_violins(axs[0, rr], subunit_exp, receptor, pcts, ds, ds_ci, kss)
    for rr, receptor in enumerate(['Kainate', 'Dopamine', '5-HT1']):
        # Prepare data
        subunit_exp, receptor, pcts, ds, ds_ci, kss = two_region_prep(subunit_data, receptor_list, receptor,
                                                                      subunit_pct_diff, subunit_d_vals,
                                                                      subunit_d_cis, subunit_ks_vals)
        # Make plot
        axs[1, rr] = make_two_violins(axs[1, rr], subunit_exp, receptor, pcts, ds, ds_ci, kss)
    fig.tight_layout(pad=3.0)
    pdf.savefig(fig)
    plt.close()
    # Third PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    for rr, receptor in enumerate(['5-HT2', '5-HT3+', 'NAalpha1']):
        # Prepare data
        subunit_exp, receptor, pcts, ds, ds_ci, kss = two_region_prep(subunit_data, receptor_list, receptor,
                                                                      subunit_pct_diff, subunit_d_vals,
                                                                      subunit_d_cis, subunit_ks_vals)
        # Make plot
        axs[0, rr] = make_two_violins(axs[0, rr], subunit_exp, receptor, pcts, ds, ds_ci, kss)
    for rr, receptor in enumerate(['NAalpha2', 'NAbeta']):
        # Prepare data
        subunit_exp, receptor, pcts, ds, ds_ci, kss = two_region_prep(subunit_data, receptor_list, receptor,
                                                                      subunit_pct_diff, subunit_d_vals,
                                                                      subunit_d_cis, subunit_ks_vals)
        # Make plot
        axs[1, rr] = make_two_violins(axs[1, rr], subunit_exp, receptor, pcts, ds, ds_ci, kss)
    fig.delaxes(axs[1][2])
    fig.tight_layout(pad=3.0)
    pdf.savefig(fig)
    plt.close()
    return pdf


def multisub_violin(subunit_data, receptor_list, pdf, receptor_median):


