
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


def adjust_country_name(country_name, adjustment='default'):
    """
    Currently very simple method to convert one country's name into that used by the GTB Report. However, likely to
    need to expand this as we work with more countries.

    Args:
        country_name: String for the original country name
        adjustment: In case multiple adjustments could be required, allow for alternative sets of name adjustments
    Returns:
        adjusted_country_name: Adjusted string
    """

    if adjustment == 'for_vaccination':
        if country_name == 'Moldova':
            return 'Republic of ' + country_name + ' (the)'
    else:
        if country_name == 'Philippines':
            return country_name + ' (the)'
        elif country_name == 'Moldova':
            return 'Republic of ' + country_name
        else:
            return country_name


def replace_specified_value(a_list, new_val, old_value):
    """
    Replace all elements of a list that are a certain value with a new value specified in the inputs.

    Args:
         a_list: The list being modified
         new_val: The value to insert into the list
         old_value: The value of the list to be replaced
    """

    return [new_val if val == old_value else val for val in a_list]




