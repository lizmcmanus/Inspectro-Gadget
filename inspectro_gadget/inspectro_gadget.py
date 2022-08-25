#!/usr/bin/env python3

"""
Tool to show the receptor expressions within an MRS region

MRS regions masks must be in MNI152 space (2 mm).

Arguments need to be given for region mask names when using command line
GABA and Glutamate regions must ine in a .tsv file "receptors.tsv"

"""

import os
import time
from inspectro_gadget.io import GadgetData, is_valid


def inspectro_gadget(mask_images, mask_labels=None, in_dir=None, out_root=None, create_overlap=False):

    # Check mask images are entered as a list
    mask_images = is_valid(mask_images, list)
    # Count number of regions
    if sum(isinstance(i, list) for i in mask_images) == 0:
        no_regions = 1
    else:
        no_regions = sum(isinstance(i, list) for i in mask_images)

    # Test if there are multiple subjects
    multi_sub = False
    if no_regions == 1:
        if len(mask_images) > 1:
            multi_sub = True
            no_subjects = len(mask_images)

    # Create mask labels if not given
    if not mask_labels:
        mask_labels = []
        for mask_no in no_regions:
            mask_labels.append(f'Region {mask_no}')

    # Create output directory
    if out_root:
        if not os.path.isdir(out_root):
            raise IsADirectoryError('The directory to create the output folder in does not exist.')
        out_dir = os.path.join(out_root, f'gadget-out_{time.strftime("%Y%m%d-%H%M%S")}')
    else:
        out_dir = os.path.join(os.getcwd(), f'gadget-out_{time.strftime("%Y%m%d-%H%M%S")}')
    if not os.path.isdir(out_dir):
        os.path.mkdir(out_dir)
    else:
        raise IsADirectoryError('Output directory already exists.')

