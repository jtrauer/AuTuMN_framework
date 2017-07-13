
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


def replace_specified_value(a_list, new_val, old_value):
    """
    Replace all elements of a list that are a certain value with a new value specified in the inputs.

    Args:
         a_list: The list being modified
         new_val: The value to insert into the list
         old_value: The value of the list to be replaced
    Return:
         List with values replaced as described
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


def prepare_denominator(list_to_prepare):
    """
    Method to safely divide a list of numbers while ignoring zero denominators by adding a very small value to any
    zeroes.

    Args:
        list_to_prepare: The list to be used as a denominator
    Returns:
        The list with zeros replaced with small numbers
    """

    return [list_to_prepare[i] if list_to_prepare[i] > 0. else 1e-10 for i in range(len(list_to_prepare))]


def elementwise_list_addition(increment, list_to_increment):
    """
    Simple method to element-wise increment a list by the values in another list of the same length
    (as is needed in output generation).

    Args:
        increment: A list of values to be added to the previous list
        list_to_increment: The original list to be incremented
    Returns:
        The resulting list with elements being the sum of the elements of the two lists
    """

    assert len(increment) == len(list_to_increment), 'Attempted to increment lists with two lists of different lengths'
    return [sum(x) for x in zip(list_to_increment, increment)]
