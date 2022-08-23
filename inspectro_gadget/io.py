#!/usr/bin/env python3

"""
Resources for handling data input, output, and organisation.
"""

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
            mrnas: list
                Which mRNA to provide information about.
        """
        self.mask_data = deepcopy(mask_data)
        self.labels = deepcopy(labels)
        self.mrnas = None
        self.no_regions = deepcopy(no_regions)
        self.no_subjects = deepcopy((no_subjects))

