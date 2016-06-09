#!/usr/bin/env python
""" This is a set of functions to calculate cis- and trans-effects. I am doing
a lot of this development in an IPython notebook. I pulled these functions out
to be used across notebooks.

"""
from __future__ import division
import pandas as pd
import numpy as np
from IPython.display import display


def marenEq(df, Eii, Eti, group):
    """ Calculates cis- and trans-effects.

    Uses the formulas from:

    Nuzhdin, S. V, Friesen, M. L., & McIntyre, L. M. (2012). Genotype-phenotype
    mapping in a post-GWAS world. Trends in Genetics: TIG, 28(9), 421-6.

    To calculate cis- and trans-effects. I initially tried to iterate over
    groups and do each calculation. This was very slow, so I vectorized the
    method. Hopefully it is still clear, though more complex. It is a lot
    faster.

    :Arguments:
        :type df: pandas.DataFrame
        :param df: DataFrame containing measures of expression for line and
            tester alleles.

        :param str Eii: The column name containing counts for the line allele.

        :param Eti: The column containing counts for the tester allele.

    :rtype: pandas.DataFrame
    :returns: A dataframe to include cis and trans effects.

    """

    # Get original column names and append the columns I want to output.
    cols = df.columns.tolist()
    cols.extend(['mu', 'cis_line', 'trans_line', 'cis_tester', 'trans_tester'])
    #cols.extend(['cis_line', 'cis_tester'])

    # Calculate the difference between tester and line in both directions. I
    # will need both of these later.
    df['diffEtiEii'] = df[Eti] - df[Eii]
    df['diffEiiEti'] = df[Eii] - df[Eti]

    # Calculate the sum of tester and line
    df['sumEtiEii'] = df[Eti] + df[Eii]

    # Group DataFrame by categories
    grp = df.groupby(group)

    # Count the number of genotypes for each exonic region
    n = grp.count().iloc[:, 0]
    n.name = 'n'

    # Calculate overall average across both alleles, for a given fusion (mu). I
    # will sum the sums and 2n.
    # mu = sum(Eti + Eii) / 2n
    mu = grp['sumEtiEii'].sum() / (2 * n)
    mu.name = 'mu'

    # Calculate sum(Eti), I need this latter to calculate Tt.
    sumEti = grp[Eti].sum()
    sumEti.name = 'sumEti'

    # Merge n, mu, and sumEti onto df
    mun = pd.concat([n, mu, sumEti], axis=1)
    mun.reset_index(inplace=True)
    munMerge = df.merge(mun, how='left', on=group)

    # Calculate (Eti - Eii / n) for Ct
    munMerge['diff_n'] = munMerge['diffEtiEii'] / munMerge['n']

    # Regroup updated df so I can sum by group again.
    grp2 = munMerge.groupby(group)

    # Tester cis-effects
    # Ct = sum(Eti - Eii / n)
    # Already calculated (Eti - Eii / n) above. Now I need to sum by group.
    Ct = grp2['diff_n'].sum()
    Ct.name = 'cis_tester'

    # Merge on Ct
    CtI = Ct.reset_index()
    ctMerge = munMerge.merge(CtI, how='left', on=group)

    # Line cis-effects
    # Ci = Eii - Eti + Ct
    ctMerge['cis_line'] = ctMerge['diffEiiEti'] + ctMerge['cis_tester']

    # Tester trans-effects
    # Tt = 2(sum(Eti)/n - mu - Ct)
    ctMerge['trans_tester'] = 2 * (ctMerge['sumEti'] / ctMerge['n'] - ctMerge['mu'] - ctMerge['cis_tester'])

    # Line trans-effects
    # Ti = 2(Eti - mu - Ct)  - Tt
    ctMerge['trans_line'] = 2 * (ctMerge[Eti] - ctMerge['mu'] - ctMerge['cis_tester']) - ctMerge['trans_tester']

    return ctMerge[cols].copy()


def marenPrintTable(df, mating_status='M', fusion_id='F10005_SI', line='sum_line', tester='sum_tester'):
    """ Pretty Maren Table. """
    # Name of columns to include in output
    out_columns = ['line', 'mating_status', 'fusion_id', 'flag_AI_combined', 'sum_both', line, tester, 'cis_line', 'cis_tester', 'trans_line', 'trans_tester', 'mean_apn']
    #out_columns = ['line', 'mating_status', 'fusion_id', 'flag_AI_combined', 'sum_both', line, tester, 'cis_line', 'cis_tester', 'mean_apn']

    # Filter and print
    out = df[(df['fusion_id'] == fusion_id) & (df['mating_status'] == mating_status)][out_columns]
    display(out.set_index(['line', 'mating_status', 'fusion_id']))
