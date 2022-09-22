#!/usr/bin/env python3

"""
Resources for plotting MRS voxel locations and the various mRNA value plots.

"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sb
from matplotlib import table
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
        cmass = np.round(center_of_mass(mask_imgs[labels[aa]]), 0).astype(int)
        # Sagittal
        mask_roi = np.ma.masked_where(mask_imgs[labels[aa]][cmass[0], :, :] == 0, mask_imgs[labels[aa]][cmass[0], :, :])
        ax = fig.add_subplot(gs[aa, 0])
        ax.imshow(bground[cmass[0], :, :].T, cmap='Greys_r', origin='lower', interpolation='none')
        ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
        ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
        ax.tick_params(left=False, bottom=False)
        # Coronal
        mask_roi = np.ma.masked_where(mask_imgs[labels[aa]][:, cmass[1], :] == 0, mask_imgs[labels[aa]][:, cmass[1], :])
        ax = fig.add_subplot(gs[aa, 1])
        ax.imshow(bground[:, cmass[1], :].T, cmap='Greys_r', origin='lower', interpolation='none')
        ax.imshow(mask_roi.T, origin='lower', interpolation='none', alpha=0.8)
        ax.set(yticklabels=[], xticklabels=[])  # remove the tick labels
        ax.tick_params(left=False, bottom=False)
        ax.set_title(labels[aa])
        # Axial
        mask_roi = np.ma.masked_where(mask_imgs[labels[aa]][:, :, cmass[2]] == 0, mask_imgs[labels[aa]][:, :, cmass[2]])
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
    PDFPages object

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
    sb.violinplot(data=receptors, x='subunit', y='values', hue='region', inner="box",
                  ax=ax, linewidth=0.1, grid_linewidth=1)
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
    the_table = table.table(ax, cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
    the_table.set_fontsize(8)
    the_table.scale(1, 0.7)
    return ax


def two_region_prep(subunit_data, receptor_list, receptor, subunit_pct_diff, subunit_d_vals, subunit_d_cis, subunit_ks_vals):
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


def two_region_radar(receptor_median, labels, receptor_list, pdf):
    """
    Create radar plots showing median gene expresion for (a) GABA and Glu subunits, and (b) neuromodulator receptor
    subunits. Shows two regions on each plot.

    Parameters
    ----------
    receptor_median: dictionary
        Dictionary with dataframes holding subunit receptor medians for each region.
    labels: list
        List with the names of each region.
    pdf: PDFPages object
        PDF object to add figures to

    Returns
    -------
    PDFPages object

    """
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2)
    # One radar for transmitters (GABA+Glu) and one for neuromodulators
    for tt, radar in enumerate(['trans_radar', 'mod_radar']):
        # Create plot
        ax = fig.add_subplot(gs[0, tt], polar=True)
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0,", "0.2", "0.4", "0.6", "0.8", "1"], color="grey", size=6)
        ax.set_ylim(0, 1)
        for rr, region in enumerate(labels):
            # Get subunit names and values
            df = receptor_median[region].loc[:, receptor_list[radar].values]
            subunits = df.columns.values
            n_subunits = len(subunits)
            values = df.iloc[0, :].values.flatten().tolist()
            values += values[:1]
            # Calculate spacing
            angles = [n / float(n_subunits) * 2 * np.pi for n in range(n_subunits)]
            angles += angles[:1]
            # Add to plot
            ax.set_xticks(angles[:-1], subunits, size=7)
            ax.set_rlabel_position(0)
            ax.plot(angles, values, linewidth=1, linestyle='solid', label=region)
            # Add legend
            # ax.legend(loc='upper left', bbox_to_anchor=(-0.4, 1), fontsize=7)
    # Add to pdf
    fig.tight_layout(pad=3.0)
    pdf.savefig(fig)
    plt.close()
    return pdf


def split_subunit_list(subunits, length):
    """
    Split a list of subunits into sublists of set length. Used to get lists of subunits that will fit on a single
    PDF page.

    Parameters
    ----------
    subunits: list
        List of subunits to split
    length: int
        How long the sublists should be

    Returns
    -------

    """
    # Loop to required length
    for i in range(0, len(subunits), length):
        # Yield means the function starts where it finished the last time
        yield subunits[i:i + length]


def multisub_prep(subunit_data, subunit, receptor_median):
    subunit_exp = []
    for subject in receptor_median.index:
        subunit_exp.append(subunit_data[subject].loc[:, subunit].dropna())
    median_dist = (receptor_median[subunit].values/np.median(receptor_median[subunit].values)*100)-100
    return subunit_exp, receptor_median[subunit], median_dist


def make_multisub_violin(ax, subunit_exp, subunit, medians, median_dist):

    subjects = medians.index.values
    n_subs = len(subjects)
    # Table contents
    rows = ['% from\nmedian']
    cell_text = [[f'{dist:0.2f}' for dist in median_dist]]
    cell_text[0].append('')  # Blank for median column
    columns = list(medians.index.values)
    columns.append('Median')
    # Plot data
    for ss in range(n_subs):
        ax.violinplot([subunit_exp[ss]], positions=[ss], showextrema=False)
        ax.scatter(n_subs, medians.iloc[ss], s=18)
    # Add overall median as dot
    ax.scatter(n_subs, np.median(medians), edgecolor='black', facecolor=None, s=24)
    ax.set_title(subunit, fontsize=6)
    ax.tick_params(axis='y', which='both', labelsize=4, width=0.5)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.grid(which='minor', linestyle='-', linewidth='0.01')
    # Add table
    the_table = table.table(ax, cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
    the_table.set_fontsize(6)
    the_table.scale(1, 1.1)
    return ax


def multisub_violin(subunit_data, receptor_list, pdf, receptor_median):
    # Split the receptor list into sublists of six to fit page
    subunit_lists = list(split_subunit_list(receptor_list.loc[receptor_list.multisub_violin, 'subunit'].values, 6))
    # Create plots
    for subunits in subunit_lists:
        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 7), linewidth=0.01)
        plt.tick_params(bottom=False, top=False, left=False, right=False)
        fig.text(0.015, 0.5, 'Normalised mRNA Expression Value', va='center', ha='center', rotation='vertical',
                 fontsize=12)
        fig.tight_layout(pad=4.0)
        for rr, subunit in enumerate(subunits[0:3]):
            # Prepare data
            subunit_exp, medians, median_dist = multisub_prep(subunit_data, subunit, receptor_median)
            # Make plot
            axs[0, rr] = make_multisub_violin(axs[0, rr], subunit_exp, subunit, medians, median_dist)
        for rr, subunit in enumerate(subunits[3:]):
            # Prepare data
            subunit_exp, medians, median_dist = multisub_prep(subunit_data, subunit, receptor_median)
            # Make plot
            axs[1, rr] = make_multisub_violin(axs[1, rr], subunit_exp, subunit, medians, median_dist)
        # Remove any empty plots
        if len(subunits) == 4:
            fig.delaxes(axs[1][1])
            fig.delaxes(axs[1][2])
        if len(subunits) == 5:
            fig.delaxes(axs[1][2])
        fig.tight_layout(pad=3.0)
        pdf.savefig(fig)
        plt.close()
    return pdf


def multisub_radar(receptor_median, receptor_list, pdf):
    """
    Create radar plot showing median gene expresion for neuromodulator receptor subunits. Shows all subjects on one
    plot.

    Parameters
    ----------
    receptor_median: dictionary
        Dictionary with dataframes holding subunit receptor medians for each region.
    receptor_list: list
        List with the names of each region.
    pdf: PDFPages object
        PDF object to add figures to

    Returns
    -------
    PDFPages object

    """
    df = receptor_median.loc[:, receptor_list.multisub_radar.values]
    subjects = df.index.values
    receptors = df.columns.values
    n_receptors = len(receptors)

    angles = [n / float(n_receptors) * 2 * np.pi for n in range(n_receptors)]
    angles += angles[:1]  # Why is this done??

    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    plt.xticks(angles[:-1], receptors, size=7)
    ax.set_rlabel_position(0)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0,", "0.2", "0.4", "0.6", "0.8", "1"], color="grey", size=6)
    plt.ylim(0, 1)

    for subject in subjects:
        values = df.loc[subject, :].values.flatten().tolist()
        values += values[:1]  # Why is this done??
        ax.plot(angles, values, linewidth=1, linestyle='solid', label=subject)
        #ax.fill(angles, values, alpha=0.1)

    # Add legend
    ax.legend(loc='upper left', bbox_to_anchor=(-0.4, 1), fontsize=7)
    fig.tight_layout(pad=3.0)
    pdf.savefig(fig)
    plt.close()
    return pdf
