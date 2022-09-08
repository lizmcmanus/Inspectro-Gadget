#!/usr/bin/env python3

"""
Resources for plotting MRS voxel locations and the various mRNA value plots.

"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sb
from scipy.ndimage import center_of_mass

def plot_masks(mask_imgs, labels, bground):
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
    return fig


def plot_overlap(overlap, labels, bground):
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
    return fig


def make_single_violin(ax, receptors, group):
    """
    Create a subplot with violin plots showing receptor expression for a group of subunits.

    Parameters
    ----------
    ax: matplotlib axis
        The axis to add the plot to
    receptors: dataframe
        Dataframe with receptor expression. Column names are subunit names
    group: string
        The subunit grouping name for use in the title

    Returns
    -------
    ax: matplotlib axis
        The axis with the plots added
    """
    sb.violinplot(data=receptors, inner="box", ax=ax, linewidth=0.1, grid_linewidth=1)
    ax.set_title(group, fontsize=6)
    #ax.set_xticks(np.arange(0, n))
    ax.tick_params(axis='x', which='major', labelsize=6, rotation=45)
    ax.tick_params(axis='y', which='major', labelsize=6)
    return ax


def single_region_violins(subunit_data, receptor_list, pdf):
    """

    Parameters
    ----------
    subunit_data
    receptor_list
    pdf

    Returns
    -------

    """
    # First PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    cnt = 0
    for rr, receptor in enumerate(['GABAA_Alpha', 'GABAA_Beta', 'GABAA_Gamma']):
        axs[0, rr] = make_single_violin(axs[0, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    for rr, receptor in enumerate(['GABAB', 'NMDA', 'AMPA']):
        axs[1, rr] = make_single_violin(axs[1, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    pdf.savefig(fig)
    # Second PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    cnt = 0
    for rr, receptor in enumerate(['mGlu(I)', 'mGlu(II)', 'mGlu(III']):
        axs[0, rr] = make_single_violin(axs[0, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    for rr, receptor in enumerate(['Kainate', 'Dopamine', '5-HT1']):
        axs[1, rr] = make_single_violin(axs[1, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    pdf.savefig(fig)
    # Third PDF page
    fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical', fontsize=12)
    fig.tight_layout(pad=4.0)
    cnt = 0
    for rr, receptor in enumerate(['5-HT2', '5-HT3+', 'NAalpha1']):
        axs[0, rr] = make_single_violin(axs[0, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    for rr, receptor in enumerate(['NAalpha2', 'NAbeta']):
        axs[1, rr] = make_single_violin(axs[1, rr],
                                        subunit_data[receptor_list.subunit[receptor_list.grouping == receptor].values],
                                        receptor)
    fig.delaxes(axs[1][2])
    pdf.savefig(fig)




