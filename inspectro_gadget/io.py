#!/usr/bin/env python3

"""
Resources for handling data input, output, and organisation.
"""

import os
import nibabel as ni
import numpy as np
import pandas as pd
from copy import deepcopy


def is_valid(var, var_type, list_type=None):
    """
    Check that the var is of a certain type.
    If type is list and list_type is specified, checks that the list contains list_type.

    Parameters
    ----------
    var: any type
        Variable to be checked.
    var_type: type
        Type the variable is assumed to be.
    list_type: type
        Like var_type, but applies to list elements.
    Returns
    -------
    var: any type
        Variable to be checked (same as input).
    Raises
    ------
    AttributeError
        If var is not of var_type
    """
    if not isinstance(var, var_type):
        raise AttributeError(f'The given variable is not a {var_type}')

    if var_type is list and list_type is not None:
        for element in var:
            _ = is_valid(element, list_type)

    return var


def test_nifti_ext(fname):
    """
    Ensure that the input filename has the expected extension for a nifti file.

    Parameters
    ----------
    fname: string
        The filename

    """
    # Check it's a valid filename
    if fname[-4:] == '.nii':
        return
    elif fname[-7:] == '.nii.gz':
        return
    else:
        raise AttributeError(f'The input file {fname} does not appear to be for a nifti image')


def load_nifti(fname):
    """
    Load data from an input nifti image.

    Parameters
    ----------
    fname: string
        Input filename

    Returns
    -------
    Image data as a numpy array

    """
    # Check it's a valid filename
    test_nifti_ext(fname)

    # Load image
    img = ni.load(fname)

    # Ensure image is correct size
    if not img.shape == (91, 109, 91):
        raise AttributeError(f'Input file {fname} does not have the correct dimensions')
    if not tuple(img.header['pixdim'][1:4]) == (2, 2, 2):
        raise AttributeError(f'Input file {fname} does not have the correct dimensions')

    return img.get_fdata()


def extract_mrna(subunit_path, region_mask):
    """
    Extract the mRNA expression values within a mask region and normalise relative to cortical expression through
    robust sigmoid method.

    Parameters
    ----------
    subunit_path: string
        Path to image file containing mRNA values
    region_mask: array
        Numpy array containing voxel mask

    Returns
    -------
    region_norm:
        Normalised expression values for all voxels within mask

    """
    mrna = ni.load(subunit_path).get_fdata()
    mrna_removed = np.array(mrna[(mrna != 0)])

    # interquartile range
    median = np.median(mrna_removed)
    Q1 = np.percentile(mrna_removed, 25)
    Q3 = np.percentile(mrna_removed, 75)
    IQRx = (Q3 - Q1) / 1.35

    # robust sigmoid
    image_norm = 1 / (1 + np.exp(-(mrna - median) / IQRx))
    region_norm = image_norm[region_mask == 1]

    return region_norm


def get_receptor_data(receptors, mask, data_dir):
    """
    Extract mRNA data for each gene from a given mask.

    mRNA expression is normalised to cortical values.

    Parameters
    ----------
    receptors: list
        List of gene names (corresponding to the names used in the image file names)
    mask: array
        Binary mask indicating the MRS voxel
    data_dir: string
        Directory where module data is located.

    Returns
    -------
    Numpy array with a column for each gene. Each row is a voxel.

    """
    out = np.zeros((len(mask[mask == 1]), len(receptors)))
    for rr, receptor in enumerate(receptors):
        out[:, rr] = extract_mrna(os.path.join(data_dir, 'mRNA_images', f'{receptor}_mirr_mRNA.nii'), mask)
    out_df = pd.DataFrame(data=out, columns=receptors)
    # Set voxels where there is no expression data to "NaN"
    out_df[out_df < 0.1] = np.nan
    return out_df


class GadgetData:
    """
    Attributes
    ----------
        mask_fnames: list
            File names for voxel masks
        labels: list
            Labels for the masks
        multi_subject: bool
            Whether the masks are from multiple subjects or not
        no_subjects: int
            How many subjects are included (defaults to 1 if it's not a multi-subject analysis)
        receptor_list: dataframe
            Dataframe containing the list of genes to be analysed and the receptor type they are related to.
        img_affine: array
            Affine matrix for the MNI template
        bground_image: array
            Array containing the background image against which masks are to be plotted.
        mask_images: dict
            Dictionary containing arrays with the region masks. Keys are region labels.
        receptor_data: dict
            Dictionary containing the normalised mRNA expression values. Keys are region labels. Arrays within have a
            row per voxel and column per gene. Multi-subject analyses contain separate arrays per subject in a list.
        ex_in_ratio: dictionary
            Excitation/inhibition ratio for each region (and participant where relevant)

    Methods
    -------

    Notes
    -----
    """

    def __init__(self, mask_fnames, labels, multi_region=False, multi_subject=False, no_subjects=1):
        """
        Initialise data object.
        """
        #  Inititate user-defined parameters
        self.mask_fnames = deepcopy(mask_fnames)
        self.labels = deepcopy(labels)
        self.no_subjects = deepcopy((no_subjects))
        self.multi_subject = deepcopy(multi_subject)
        self.multi_region = deepcopy(multi_region)

        # Load built-in data
        data_dir = '/media/storage2/Inspectro-Gadget/inspectro_gadget/data' ###### Change this!!!
        #data_dir = os.path.dirname(inspect.getfile(inspectro_gadget))
        self.receptor_list = pd.read_csv(os.path.join(data_dir, 'GroupedReceptors.tsv'), delimiter='\t', header=0)
        tmp_img = ni.load(os.path.join(data_dir, 'MNI152_T1_2mm.nii.gz'))
        self.img_affine = tmp_img.affine
        self.bground_image = tmp_img.get_fdata()

        # Load mask(s) and receptor data
        self.mask_images = {}
        self.receptor_data = {}
        if self.multi_subject:
            for ss in range(no_subjects):
                print(f'Loading data for subject {self.labels[ss]}')
                self.mask_images[self.labels[ss]] = load_nifti(self.mask_fnames[ss])
                self.receptor_data[labels[ss]] = get_receptor_data(self.receptor_list.iloc[:, 0].values,
                                                                       self.mask_images[labels[ss]], data_dir)
        elif self.multi_region:
            for ll, label in enumerate(labels):
                self.mask_images[label] = load_nifti(self.mask_fnames[ll][0])
                self.receptor_data[label] = get_receptor_data(self.receptor_list.iloc[:, 0].values,
                                                              self.mask_images[label], data_dir)
        else:
            self.mask_images[labels[0]] = load_nifti(self.mask_fnames[0])
            self.receptor_data[labels[0]] = get_receptor_data(self.receptor_list.iloc[:, 0].values,
                                                          self.mask_images[labels[0]], data_dir)


        # Inititate common stats variables
        self.ex_in_ratio = {}
        self.receptor_median = {}
        self.overlap_image = {}
