#!/usr/bin/env python3

"""
Tool to show the receptor expressions within an MRS region

MRS regions masks must be in MNI152 space (2 mm).

Arguments need to be given for region mask names when using command line
GABA and Glutamate regions must ine in a .tsv file "receptors.tsv"

"""

import os
import time
from inspectro_gadget.io import GadgetData, load_nifti
from inspectro_gadget import plotting, stats
from matplotlib.backends.backend_pdf import PdfPages


def gadget(mask_fnames, mask_labels=None, out_root=None, create_overlap=False, bground_fname=None):
    """
    Main function that runs the analysis. User provides the location(s) of the masks to be used and the relevant output
    will be produced depending upon whether masks for one region, two regions, or multiple subjects are entered.

    Parameters
    ----------
    mask_fnames: list
        List containing mask filename(s). If one region then a list with a single string. If two regions then a list
        with two filenames, each in their own sub-list. If multiple subjects then a list with multiple filenames.
    mask_labels: list
        List containing the labels for each mask. Structure follows that used for filenames.
    out_root: string
        Directory to which output should be written. If none is given then a new directory will be made in the current
        working directory.
    create_overlap: Boolian
        Whether to create a nifti file with the subject overlap in a multi-subject analysis.
    bground_fname: string
        Filename of alternative background image for mask location plots. Must be in MNI152 2mm space.

    Returns
    -------

    """

    # Check mask images are entered as a list
    mask_fnames = is_valid(mask_fnames, list)
    # Count number of regions
    if sum(isinstance(i, list) for i in mask_fnames) == 0:
        multi_region = False
    else:
        multi_region = True
        no_regions = sum(isinstance(i, list) for i in mask_fnames)
        no_subjects = 1

    # Test if there are multiple subjects
    multi_sub = False
    if not multi_region:
        if len(mask_fnames) > 1:
            multi_sub = True
            no_subjects = len(mask_fnames)

    # Create mask labels if not given
    if not mask_labels:
        mask_labels = []
        for mask_no in range(1, no_regions+1):
            mask_labels.append(f'Region {mask_no}')

    # Create output directory
    if out_root:
        if not os.path.isdir(out_root):
            raise IsADirectoryError('The directory to create the output folder in does not exist.')
        out_dir = os.path.join(out_root, f'gadget-out_{time.strftime("%Y%m%d-%H%M%S")}')
    else:
        out_dir = os.path.join(os.getcwd(), f'gadget-out_{time.strftime("%Y%m%d-%H%M%S")}')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    else:
        raise IsADirectoryError('Output directory already exists.')

    # Collect all required data
    data = GadgetData(mask_fnames, mask_labels, multi_region=multi_region, multi_subject=multi_sub,
                      no_subjects=no_subjects)

    # Replace background image if the user provides one
    if bground_fname:
        data.bground_image = load_nifti(bground_fname)

    # Calculate statistics
    if data.multi_subject:
        data.receptor_median = stats.subject_median(data.receptor_data, data.receptor_list)
        for subject in data.labels:
            data.ex_in_ratio[subject] = stats.ex_in(data.receptor_data[subject], data.receptor_list)
    else:
        for region in data.labels:
            data.receptor_median[region] = stats.region_median(data.receptor_data[region], data.receptor_list)
            data.ex_in_ratio[region] = stats.ex_in(data.receptor_data[region], data.receptor_list)
        if data.multi_region:
            data.subunit_d_vals, data.subunit_d_cis, data.subunit_pct_diff, data.subunit_ks_vals = stats.compare_regions(data.receptor_data,
                                                                                                                         data.receptor_list)
    # Create output PDF
    with PdfPages('gadget-output.pdf') as pdf:
        # Mask images
        if data.multi_subject:
            pdf = plotting.plot_masks(data.overlap_image, data.labels, data.bground_image, pdf)
        else:
            pdf = plotting.plot_masks(data.mask_images, data.labels, data.bground_image, pdf)

        # Violin plots
        if data.multi_subject:
            pdf = plotting.multi_sub_violin(data.receptor_data, data.receptor_list, pdf, data.receptor_median)
        else:
            if data.multi_region:
                pdf = plotting.two_region_violins(data.receptor_data, data.receptor_list, pdf, data.subunit_pct_diff,
                                                  data.subunit_d_vals, data.subunit_d_cis, data.subunit_ks_vals)
            else:
                pdf = plotting.single_region_violins(data.receptor_data[data.labels[0]], data.receptor_list, pdf)

        # Radar plots
        if data.multi_subject:
            pdf = plotting.multisub_radar(data.receptor_median, data.receptor_list, pdf)
        if data.multi_region:
            pdf = plotting.two_region_radar(data.receptor_median, data.labels, pdf)
