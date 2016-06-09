#!/usr/bin/env python

import seaborn as sns
import pandas as pd
from sklearn.neighbors import DistanceMetric
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def buildPalette(df):

    levels = df['group'].value_counts().index.tolist()

    clist = sns.color_palette("Set2", len(levels))

    colors = dict()
    p = list()

    for i, lvl in enumerate(levels):
        p.append(mpatches.Patch(color=clist[i], label=lvl))

        for geno in df.loc[df['group'] == lvl, 'line'].tolist():
            colors[geno] = clist[i]

    return colors, p


def cluster(ai, maren, fusion):
    """ """

    # Filter out current fusion of interest
    ai = ai[ai['fusion_id'] == fusion]

    maren = maren[maren['fusion_id'] == fusion]
    maren.set_index('fusion_id', inplace=True)

    # Calculate Euclidean Distance of cis-line effects
    euc = DistanceMetric.get_metric('euclidean')
    dist = pd.DataFrame(euc.pairwise(maren.T), columns=maren.columns, index=maren.columns)

    # Associate colors with different group (line, tester, noai, ambiguous)
    colors, patches = buildPalette(ai)
    c = [colors[x] for x in maren.columns]

    # Draw plot
    fig = sns.clustermap(dist, row_colors=c)
    fig.ax_heatmap.legend(handles=patches, loc=(1, 1.1))
    fig.ax_col_dendrogram.set_title(fusion, fontsize=18)

    return fig
