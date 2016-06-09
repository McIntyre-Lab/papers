#!/usr/bin/env python
""" Set of functions for plotting things in the ASE project """
import pandas as pd
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import inspect


def pdPlotDefaults():
    """ Returns a dictionary of all of the pandas plot default values. """

    # Lookat pandas plot method and get options
    options = inspect.getargspec(pd.DataFrame.plot)

    # Return only kward options with their defaults
    return dict(zip(options.args[-len(options.defaults):], options.defaults))


def dfPanelScatter(df, x, y, group, figsize=(20, 20), colorCol=False, colormap='rainbow', vline=False, hline=False, plot_title=False, **kwargs):
    """ Plot a panel of scatter plots from a pandas DataFrame.

    :Arguments:

        :type df: pandas.DataFrame
        :param df: Dataframe to plot.

        :param str x: Column name to plot on x-axis.

        :param str y: Column name to plot on y-axis.

        :param str|list group: Column name with group information.

        :param str colorCol: Column name with numeric values for coloring.

        :param str colormap: Name of the colormap to use if plotting colors.

        :param dict **kwargs: Additional pandas.DataFrame.plot keyword
            arguments.

    :rtype: mp.pyplot.Figure
    :returns: Generates a plot and returns the

    """
    # Get pandas plot defaults
    opts = pdPlotDefaults()

    # Update Defaults
    opts.update({'x': x, 'y': y, 'kind': 'scatter'})
    if colorCol:
        opts.update({'c': colorCol, 'colormap': colormap})

    opts.update(kwargs)

    # Group and iterate
    grp = df.groupby(group)

    # How many plots are there and what is the dimensions of the panel
    n = len(grp.groups.keys())
    dim = int(np.ceil(np.sqrt(n)))
    fig, axes = plt.subplots(dim, dim, figsize=(20, 20))
    axes = axes.ravel()

    if plot_title:
        fig.suptitle(plot_title, fontsize=18)

    cnt = 0
    for i, dfG in grp:
        ax = axes[cnt]
        opts.update({'title': i, 'ax': ax})
        dfG.plot(**opts)
        if vline:
            ax.axvline(vline)
        if hline:
            ax.axhline(hline)
        cnt += 1


# Function to plot theta vs cis-line
def plotComparison(df, theta='q5_mean_theta', ci='cis_line', seed=False):
    """ Function to plot theta's vs Ci's.

    I am trying to figure out if the measures of theta is concordant to the
    Ci's.  This function generates plots to look at concordance.

    :type df: pandas.DataFrame
    :param str theta: name of the column with thetas
    :param stt Ci: name of the column with cis-effect for the line
    """
    if seed:
        np.random.seed(seed)

    dfI = df.reset_index()
    grp = dfI.groupby('fusion_id')
    fus = grp.groups.keys()
    randomFus = np.random.choice(fus, 25)

    fig, axes = plt.subplots(5, 5, figsize=(20, 20))
    axes = axes.ravel()

    for i, val in enumerate(randomFus):
        dfI[dfI['fusion_id'] == val].plot(kind='scatter', x=ci, y=theta, ax=axes[i], title=val, rot=90, c='flag_AI_combined')
        axes[i].axhline(0.5, color='r', lw=2)
        axes[i].set_ylim(0, 1)
    plt.tight_layout()
