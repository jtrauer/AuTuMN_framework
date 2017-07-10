
"""
Note that this module is intended only to contain stand-alone functions for use by multiple other modules.
Object-oriented structures are not intended to be kept here.
"""


def find_first_list_element_at_least_value(list, value):
    """
    Simple method to return the index of the first element of a list that is greater than a specified value.

    Args:
        list: The list of floats
        value: The value that the element must be greater than
    """

    return next(x[0] for x in enumerate(list) if x[1] >= value)

