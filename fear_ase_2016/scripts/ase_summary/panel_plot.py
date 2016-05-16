#!/usr/bin/env python
import matplotlib.pyplot as plt


def plotPanelKDE(ms, dat):
    """ Function to plot a panel of distributions of $\theta$. """
    # mating status
    if ms == 'M':
        MS = 'Mated'
    else:
        MS = 'Virgin'

    # Group by line
    grp = dat.groupby('line')

    # Plot all distributions and color bg by APN rank
    fig, axes = plt.subplots(8, 9, figsize=(20, 20))
    fig.suptitle(MS, fontsize=18)
    axes = axes.ravel()

    # Iterate over lines and plot each distribution.
    cnt = 0
    for i, val in grp:
        ax = axes[cnt]
        val[['q4_mean_theta', 'q5_mean_theta', 'q6_mean_theta']].plot(kind='kde', ax=ax, color=['b', 'r', 'g'], legend=False, title=i)
        ax.axvline(0.5, color='k')
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        cnt += 1
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])


def plotPanelSnpIndel(ms, dat):
    """ Function to plot a panel of SNP/INDEL counts vs $\theta$. """
    # mating status
    if ms == 'M':
        MS = 'Mated'
    else:
        MS = 'Virgin'

    # Group by line
    grp = dat.groupby('line')

    # Plot all distributions and color bg by APN rank
    fig, axes = plt.subplots(8, 9, figsize=(20, 20), sharey=True, sharex=True)
    fig.suptitle(MS, fontsize=18)
    axes = axes.ravel()

    # Iterate over lines and plot each distribution.
    cnt = 0
    for i, val in grp:
        ax = axes[cnt]
        val.plot(kind='scatter', x='num_snps', y='q5_mean_theta', ax=ax, color='b', label='SNPs', legend=False, title=i)
        val.plot(kind='scatter', x='num_indels', y='q5_mean_theta', ax=ax, color='r', marker='^', label='INDELs', legend=False)
        ax.axhline(0.5, color='r')
        ax.set_ylabel(r'$\theta$')
        ax.set_xlabel('Number Polymorphisms')
        ax.get_xaxis().set_ticks([])
        ax.set_ylim(-0.2, 1.2)
        ax.set_xlim(-0.2, 100)
        cnt += 1
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])


def plotPanelAPN(ms, dat):
    """ Function to plot a panel of SNP/INDEL counts vs $\theta$. """
    # mating status
    if ms == 'M':
        MS = 'Mated'
    else:
        MS = 'Virgin'

    # Group by line
    grp = dat.groupby('line')

    # Plot all distributions and color bg by APN rank
    fig, axes = plt.subplots(8, 9, figsize=(20, 20), sharey=True, sharex=True)
    fig.suptitle(MS, fontsize=18)
    axes = axes.ravel()

    # Iterate over lines and plot each distribution.
    cnt = 0
    for i, val in grp:
        ax = axes[cnt]
        val.plot(kind='scatter', x='mean_apn', y='q5_mean_theta', ax=ax, color='b', legend=False, title=i, rot=45)
        ax.axhline(0.5, color='r')
        ax.set_ylabel(r'$\theta$')
        ax.set_xlabel('Mean APN')
        ax.set_ylim(-0.2, 1.2)
        ax.set_xlim(0, 10000)
        cnt += 1
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
