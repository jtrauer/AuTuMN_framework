
from __future__ import print_function
from xlrd import open_workbook
from numpy import nan
import numpy
import os
import tool_kit


#######################################
###  Individual spreadsheet readers ###
#######################################

class GlobalTbReportReader:
    """
    Reader object for the WHO's Global TB Report 2016. Illustrates general structure for spreadsheet readers.
    """

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
        self.year_indices = {}
        self.country_to_read = tool_kit.adjust_country_name(country_to_read)

    def parse_col(self, col):
        """
        Read and interpret a column of the spreadsheet

        Args:
            col: The column to be read
        """

        col = tool_kit.replace_specified_value(col, nan, '')

        # if it's the country column (the first one), find the indices for the country being simulated
        if col[0] == 'country':
            for i in range(len(col)):
                if col[i] == self.country_to_read: self.indices += [i]

        # ignore irrelevant columns
        elif 'iso' in col[0] or 'g_who' in col[0] or 'source' in col[0]:
            pass

        # find years to read from year column
        elif col[0] == 'year':
            for i in self.indices:
                self.year_indices[int(col[i])] = i

        # get data from the remaining (data) columns
        else:
            self.data[str(col[0])] = {}
            for year in self.year_indices:
                if not numpy.isnan(col[self.year_indices[year]]):
                    self.data[col[0]][year] = col[self.year_indices[year]]

    def get_data(self):
        """
        Return the read data.
        """

        return self.data


#########################
###  Master functions ###
#########################

def read_xls_with_sheet_readers(sheet_readers):
    """
    Runs each of the individual readers (currently only one) to gather all the data from the input spreadsheets.

    Args:
        sheet_readers: The sheet readers that have been collated into a list
    Returns:
        All the data from the reading process as a single object
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

        # if the workbook was found to be available available, read the sheet in question
        else:
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
        sheets_to_read: A list containing the strings that are also the 'keys' attribute of each reader
        country: Country being read

    Returns:
        A single data structure containing all the data to be read
    """

    sheet_readers = []
    if 'tb' in sheets_to_read: sheet_readers.append(GlobalTbReportReader(country))
    for reader in sheet_readers: reader.filename = os.path.join(reader.filename)
    return read_xls_with_sheet_readers(sheet_readers)


