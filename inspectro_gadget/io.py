#!/usr/bin/env python3

"""
Resources for handling data input, output, and organisation.
"""

import nibabel as ni
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




class GadgetData():
    """
    Attributes
    ----------

    Methods
    -------

    Notes
    -----
    """

    def __init__(self, mask_data, labels, mrnas, no_regions, no_subjects):
        """
        Initialise data object.

            Parameters
            ----------
            mask_data: list
                A list comtaining numpy arrays of the MRS masks to work with
                If a single mask then a list with a single 3d array;
                if multiple masks a list of multiple 3d;
                if subject masks a list of many strings
            labels: list
                Labels for the masks


        """
        self.mask_data = deepcopy(mask_data)
        self.labels = deepcopy(labels)
        self.mrnas = None
        self.no_regions = deepcopy(no_regions)
        self.no_subjects = deepcopy((no_subjects))

