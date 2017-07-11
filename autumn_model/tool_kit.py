
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


def sort_dictionary_to_list_by_keys(dictionary):
    """
    Extract a list of paired tuples with each tuple being a key, value pair from the input dictionary and the order of
    the list according to the order of the keys of the dictionary.

    Args:
         dictionary: Input dictionary with keys that have an inherent ordering
    Return:
         Ordered list containing pairs that were keys and values
    """

    return [(i, dictionary[i]) for i in sorted(dictionary)]


def extract_list_of_paired_tuples_to_list(list_of_paired_tuples):
    """
    Extract the elements of the list created by sort_dictionary_to_list_by_keys into two similarly ordered lists.

    Args:
         list_of_paired_tuples: Expected to be the output of sort_dictionary_to_list_by_keys
    Return:
         Two lists, one for the keys of the original dictionary and the other for the values
    """

    keys = [i[0] for i in list_of_paired_tuples]
    values = [i[1] for i in list_of_paired_tuples]
    return keys, values


def extract_ordered_list_from_dictionary(dictionary):
    """
    Simple function to just combine the two preceding functions.

    Return:
          Return the two ordered lists directly
    """

    return extract_list_of_paired_tuples_to_list(sort_dictionary_to_list_by_keys(dictionary))



