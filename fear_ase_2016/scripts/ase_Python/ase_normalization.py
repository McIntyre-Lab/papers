#!/usr/bin/env python
""" These are different functions for normalizing data. """
from __future__ import division
import pandas as pd
import numpy as np


def meanCenter(df, columns, group):
    """Mean centering.

    Mean center allele specific counts by using both line and tester counts.

    :Arguments:
        :type df: pandas.DataFrame
        :param df: DataFrame containing the value and group.

        :param list columns: Columns with allele specific counts

        :param str group: The name of the column with group information.

    :rtype: pandas.DataFrame
    :returns: Updates original DataFrame with normalized values added.

    """

    # I am mean centering allele counts. To do this I want to take the mean of
    # all the allele counts (i.e., line and tester together). Then center each
    # separately.
    grp = df.groupby(group)
    centered = grp[columns].apply(lambda x: x - np.mean(x.values))
    centered.columns = ['mean_center_' + x for x in columns]

    df = df.join(centered)
    return df


def meanStd(df, column, group):
    """ Function to calculate mean standardization.

    :Arguments:
        :type df: pandas.DataFrame
        :param df: DataFrame containing the value and group to standardize.

        :param str column: The name of the column with values you want standardize.

        :param str group: The name of the column with group information.

    :rtype: pandas.DataFrame
    :returns: Updates original DataFrame with normalized values added.

    """
    # Formula for mean standardization.
    meanStandardize = lambda x: (x - x.mean()) / x.std()

    grp = df.groupby(group)
    df['mean_std_' + column] = grp[column].transform(meanStandardize)


def q3Norm(sf):
    """ Calculate the upper quartile normalization.

    :Arguments:
        :type sf: pandas.Series
        :param sf: Series with the column of interes

    :rtype: pandas.Series
    :returns: Series with the upper quartile normalized value.

    """
    # Calculate upper quantile for each mating_status*genotype
    q3 = sf.groupby(level=[0, 1]).quantile(q=0.75)

    # Calculate median q3 by mating_status
    medQ3 = q3.groupby(level=0).median()

    # Function to calculate normalization
    norm = lambda x: x / q3 * medQ3

    # Normalize and return
    # I need to unstack the fusions into wide for calculation and then stack them
    # back.
    return sf.unstack(level=-1).apply(norm).stack(dropna=True)
