""" A class for generating flags.

Tables of binary flags are commonly used for looking at frequencies and filter
rows. Here is a class to quickly generate these types of tables.

"""
# Imports
import numpy as np
import pandas as pd
from itertools import combinations
from itertools import product
import scipy.stats as stats
from IPython.display import display, HTML


def in_ipynb():
    """ Function to test if script is running in IPython. """
    try:
        if get_ipython().config:
            return True
    except NameError:
        return False


def pprint(obj, tag='b'):
    """ Pretty print an object.

    If running in an ipython, use HMTL to make nicer output.
        * If DataFrame, make HTML table.
        * If string then surround with html tag, default is bold.

    :Arguments:
        :param obj: Any type of object that is affected by IPythons HTML
            display.

        :param str tag: An HTML tag to surround strings with, default is 'b'
            for bold.

    """
    if in_ipynb():
        if isinstance(obj, str):
            display(HTML('<{}>{}</{}>'.format(tag, obj, tag)))
        else:
            display(obj)
    else:
        print(obj)


class FlagsDataFrame(pd.DataFrame):
    """ Create subclass of pandas.DataFrame to add some functionality. """
    def __init__(self, *args, **kwargs):
        """ Creates a FlagsDataFrame.

        Optional pd.DataFrame *args and **kwargs can be passed, but no
        arguments are required.

        :returns: A FlagsDataFrame object
        :rtype: FlagsDataFrame
        """
        super(FlagsDataFrame, self).__init__(*args, **kwargs)

    def addColumn(self, column, mask=[]):
        """ Add a column and set values to 0.

        Adds a column to the flag data frame and sets the values to 0. If mask
        is present it sets True values in the mask to 1 in the data frame.

        :Arguments:
            :param str column: Name of the column to add.

            :param mask: A boolean mask.
            :type mask: pd.Series

        :returns: Updates DataFrame inplace with new column.

        """
        self[column] = 0

        if len(mask) > 0:
            self.updateMask(mask=mask, column=column)

    def updateMask(self, mask, column=''):
        """ Update the dataframe with 1's if the mask value is true.

        :Arguments:
            :param mask: List of mask values.
            :type mask: list of pd.Series of bool

            :param str column: Column name to update in the flag frame. Not required

        :returns: Updated FlagsDataFrame inplace.

        """
        # Update the values to 1's if they are true in the mask
        if len(column) > 0:
            self.loc[mask.index, column] = mask.astype(int)
        else:
            self.loc[mask.index, self.columns] = mask.astype(int)

    def _runCrosstab(self, index, columns, margins=True, *args, **kwargs):
        """ Uses Pandas crosstab to generate crosstabulation tables.

        Creates a cross tabulation table with pandas.crosstab and calculates
        fisher's exact test and a chi^2 test using scipy.

        :Arguments:
            :param str index: Name of the first column to pass to pd.crosstab

            :param str columns: Name of the second column to pass to pd.crosstab

            :param margins: Option to output margins.

            :param *args: other options to pass to pd.crosstab
            :param **kwargs: other named options to pass to pd.crosstab

        """
        ct = pd.crosstab(index, columns, margins=margins, *args, **kwargs)
        ct.rename(index={'All': 'Total'}, columns={'All': 'Total'}, inplace=True)

        # Calculate statistics
        ## If margins included then remove for calculation
        if margins:
            twoByTwo = ct.iloc[:-1, :-1]
        else:
            twoByTwo = ct

        try:
            ## Get fisher's and chi square
            fisher = stats.fisher_exact(twoByTwo)
            chi = stats.chi2_contingency(twoByTwo, correction=False)
            tests = pd.DataFrame(data=[fisher, chi[0:2]], index=['Fisher-Exact (two-sided)', 'Chi-square'], columns=['value', 'Prop'])

            ## Display
            pprint('{} vs {}'.format(index.name, columns.name), tag='h3')
            pprint(ct)
            pprint(tests)
            print('')    # print blank line
        except:
            pprint('{} vs {}'.format(index.name, columns.name), tag='h3')
            pprint(ct)
            print('Could not calculate tests, check your table')
            print('')   # print blank line

    def flagCrosstab(self, columns=None, *args, **kwargs):
        """ Calculate crosstab and stats on FlagsDataFrame.

        Calculates the crosstab and enrichment stats on the current DataFrame
        and prints results.

        :Arguments:
            :param columns: Name of columns to calculate crosstab. If no
                columns are provided then all pairwise combinations will be
                done. If a single column is provided then it will be compared
                to all other columns. If a list of columns is provided only
                those columns will be compared.

            :param *args: Additional arguments to pass to pd.crosstab
            :param **kwargs: Additional arguments to pass to pd.crosstab
        """
        if columns is None:
            # If not column names are given, then do all pairwise combinations
            for combo in combinations(self.columns, 2):
                flag1 = self[combo[0]]
                flag2 = self[combo[1]]
                self._runCrosstab(flag1, flag2, *args, **kwargs)

        elif isinstance(columns, str):
            # If 1 column name given as a string, then compare it to all other
            # columns
            for combo in product([columns], self.columns):
                if combo[0] != combo[1]:
                    flag1 = self[combo[0]]
                    flag2 = self[combo[1]]
                    self._runCrosstab(flag1, flag2, *args, **kwargs)

        elif hasattr(columns, '__iter__'):
            # If a list or array is given then do the following
            if len(columns) == 1:
                # If 1 column given in a list then treat it like if a sting was
                # given
                for combo in product(columns, self.columns):
                    if combo[0] != combo[1]:
                        flag1 = self[combo[0]]
                        flag2 = self[combo[1]]
                        self._runCrosstab(flag1, flag2, *args, **kwargs)

            elif len(columns) == 2:
                # If two column names provided only do a cross tab for those
                flag1 = self[columns[0]]
                flag2 = self[columns[1]]
                self._runCrosstab(flag1, flag2, *args, **kwargs)

            elif len(columns) > 2:
                # If more than two columns given then do those pairwise combos
                for combo in combinations(columns, 2):
                    flag1 = self[combo[0]]
                    flag2 = self[combo[1]]
                    self._runCrosstab(flag1, flag2, *args, **kwargs)

    @property
    def _constructor(self):
        return FlagsDataFrame

    @property
    def _constructor_sliced(self):
        return FlagsSeries


class FlagsSeries(pd.Series):
    @property
    def _constructor(self):
        return FlagsSeries

    @property
    def _constructor_expanddim(self):
        return FlagsDataFrame
