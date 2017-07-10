# -*- coding: utf-8 -*-

from __future__ import print_function
from xlrd import open_workbook
from numpy import nan
import numpy
import os
import tool_kit


def is_all_same_value(a_list, test_val):
    """
    Simple method to find whether all values in list are equal to a particular value.

    Args:
        a_list: The list being interrogated
        test_val: The value to compare the elements of the list against
    """

    for val in a_list:
        if val != test_val: return False
    return True


def replace_specified_value(a_list, new_val, old_value):
    """
    Replace all elements of a list that are a certain value with a new value specified in the inputs.

    Args:
         a_list: The list being modified
         new_val: The value to insert into the list
         old_value: The value of the list to be replaced
    """

    return [new_val if val == old_value else val for val in a_list]


def parse_year_data(year_data, blank, endcolumn):
    """
    Code to parse rows of data that are years.

    Args:
        year_data: The row to parse
        blank: A value for blanks to be ignored
        endcolumn: Column to end at
    """

    year_data = replace_specified_value(year_data, nan, blank)
    assumption_val = year_data[-1]
    year_vals = year_data[:endcolumn]
    if is_all_same_value(year_vals, nan):
        return [assumption_val] 
    else:
        return year_vals


#######################################
###  Individual spreadsheet readers ###
#######################################


class GlobalTbReportReader:

    def __init__(self, country_to_read):

        self.data = {}
        self.tab_name = 'TB_burden_countries_2016-04-19'
        self.key = 'tb'
        self.parlist = []
        self.filename = 'xls/gtb_data.xlsx'
        self.start_row = 1
        self.horizontal = False
        self.start_column = 0
        self.indices = []
        self.country_to_read = tool_kit.adjust_country_name(country_to_read)

    def parse_col(self, col):

        col = replace_specified_value(col, nan, '')

        # if it's the country column (the first one)
        if col[0] == 'country':

            # find the indices for the country in question
            for i in range(len(col)):
                if col[i] == self.country_to_read: self.indices += [i]

        elif 'iso' in col[0] or 'g_who' in col[0] or 'source' in col[0]:
            pass

        elif col[0] == 'year':
            self.year_indices = {}
            for i in self.indices:
                self.year_indices[int(col[i])] = i

        # all other columns
        else:
            self.data[str(col[0])] = {}
            for year in self.year_indices:
                if not numpy.isnan(col[self.year_indices[year]]):
                    self.data[col[0]][year] = col[self.year_indices[year]]

    def get_data(self):

        return self.data


######################
### Master scripts ###
######################


def read_xls_with_sheet_readers(sheet_readers):
    """
    Runs the individual readers to gather all the data from the sheets

    Args:
        sheet_readers: The sheet readers that were previously collated into a list
    Returns:
        All the data for reading as a single object
    """

    result = {}
    for reader in sheet_readers:

        # check that the spreadsheet to be read exists
        try:
            print('Reading file', os.getcwd(), reader.filename)
            workbook = open_workbook(reader.filename)

        # if sheet unavailable, print error message but continue
        except:
            print('Unable to open spreadsheet')

        else:
            # if the workbook is available, read the sheet in question
            sheet = workbook.sheet_by_name(reader.tab_name)

            # read in the direction that the reader expects (either horizontal or vertical)
            if reader.horizontal:
                for i_row in range(reader.start_row, sheet.nrows):
                    reader.parse_row(sheet.row_values(i_row))
            else:
                for i_col in range(reader.start_column, sheet.ncols):
                    reader.parse_col(sheet.col_values(i_col))
            result[reader.key] = reader.get_data()

    return result


def read_input_data_xls(sheets_to_read, country=None):
    """
    Compile sheet readers into a list according to which ones have been selected.
    Note that most readers now take the country in question as an input,
    while only the fixed parameters sheet reader does not.

    Args:
        from_test: Whether being called from the directory above
        sheets_to_read: A list containing the strings that are also the
            'keys' attribute of the reader
        country: Country being read for

    Returns:
        A single data structure containing all the data to be read
            (by calling the read_xls_with_sheet_readers method)
    """

    sheet_readers = []
    if 'tb' in sheets_to_read:
        sheet_readers.append(GlobalTbReportReader(country))
    for reader in sheet_readers:
        reader.filename = os.path.join(reader.filename)
    return read_xls_with_sheet_readers(sheet_readers)


