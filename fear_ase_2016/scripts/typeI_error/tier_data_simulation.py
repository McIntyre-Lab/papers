#!/usr/bin/env python

# Built-in packages
import os
import argparse

# Add-on packages
import numpy as np
from scipy.stats import poisson
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def getOptions():
    """ Function to pull in arguments """
    description = """ Script to run Type I Error simulation"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", dest="fname", action='store', required=True, help="Name of raw counts file. [Required]")
    parser.add_argument("--output", dest="oname", action='store', required=True, help="Name of simulated output file. [Required]")
    parser.add_argument("--fig", dest="figname", action='store', required=True, help="PDF file to store figures. [Required]")

    parser.add_argument("--number-genes", dest="nGenes", action='store', default=10 ** 4, type=int, required=False, help="Number of genes to simulate [Default: 10^4]")

    parser.add_argument("--bias", dest="bias", action='store_true', required=False, help="If given then model will be run with bias [Optional].")
    parser.add_argument("--AI", dest="ai", action='store_true', required=False, help="If given then model will be run with AI [Optional].")

    parser.add_argument("--qtrue", dest="qTrue", action='store', default=0.45, type=float, required=False, help="Set the value of qTrue [Default: 0.45]. Note if --bias is not set then qTrue is 0.5")
    parser.add_argument("--half", dest="half", action='store_true', required=False, help="Set 1/2 of simulated genes to qtrue and the other half to 1-qtrue [Default]")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.")

    args = parser.parse_args()

    # If bias option not given then set qTrue to 0.5
    if not args.bias:
        args.qTrue = 0.5

    return args


def runDNASim(out, dat, dna, sampleIDs):
    """ """
    sum1 = dat[dna[0]].mean(axis=1)
    sum2 = dat[dna[1]].mean(axis=1)
    sumMean = (sum1 + sum2)

    qTrue = 0.5
    R = 1
    b = 1 / ((R - 1) * qTrue + 1)
    a = R * b

    # Simulate DNA Reads for Sample 1
    out['DNATrueMean1'] = a * sumMean
    lambda1 = (1 - qTrue) * out['DNATrueMean1']
    random1 = lambda1.apply(lambda x: pd.Series(poisson.rvs(x, size=3)))
    dna1 = [x.replace('RNA', 'DNA') for x in sampleIDs[0]]
    out[dna1] = random1

    # Simulate DNA Reads for Sample 2
    out['DNATrueMean2'] = b * sumMean
    lambda2 = qTrue * out['DNATrueMean2']
    random2 = lambda2.apply(lambda x: pd.Series(poisson.rvs(x, size=3)))
    dna2 = [x.replace('RNA', 'DNA') for x in sampleIDs[1]]
    out[dna2] = random2

    # Create diagnostic plots
    return checkPlots(out, (dna1, dna2), sumMean, R, 'DNA')


def runRNASim(out, sumMean, sampleIDs, R):
    """ Simulate reads from a Poisson distribution.

    Get sample1 AND sample2 true means where:

        qTrue * a + (1 - qTrue) * b = 1
        a / b = R

    Args:
        out: pd.DataFrame that contains information I want to output.

        sumMean: pd.Series that contains the sum of the average number of reads
            from Sample1 and Sample2.

        SampleIDs: tuple of lits that contain the sample name column headers.

        R: is the ratio of mean1 / mean2 reads, when R = 1 then no AI, when R > 1
            then AI towards sample1 after correcting for bias.

    Updates:
        out: pd.DataFrame that contains information I want to output.

    """
    b = 1 / ((R - 1) * (1 - out['qTrue']) + 1)
    a = R * b

    # Simulate Reads for Sample 1
    out['trueMean1'] = a * sumMean
    lambda1 = (1 - out['qTrue']) * out['trueMean1']
    random1 = lambda1.apply(lambda x: pd.Series(poisson.rvs(x, size=3)))
    out[sampleIDs[0]] = random1

    # Simulate Reads for Sample 2
    out['trueMean2'] = b * sumMean
    lambda2 = out['qTrue'] * out['trueMean2']
    random2 = lambda2.apply(lambda x: pd.Series(poisson.rvs(x, size=3)))
    out[sampleIDs[1]] = random2

    # Create diagnostic plots
    return checkPlots(out, sampleIDs, sumMean, R, 'RNA')


def checkPlots(out, sampleIDs, sumMean, R, title):
    """ Create a couple of diagnostic plots to see how well simulations did. """
    sum1 = out[sampleIDs[0]].mean(axis=1)
    sum2 = out[sampleIDs[1]].mean(axis=1)

    mask = (sum1 > 0) & (sum2 > 0)
    ot = out[mask]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle(title)
    ax1.plot((sum1[mask] / (1 - ot['qTrue'])) / (sum2[mask] / ot['qTrue']))
    ax1.set_ylabel('Line / Tester after bias correction')
    ax1.set_xlabel('Simulated Exonic Regions')
    ax1.axhline(R, color='r')

    x = np.log10(sumMean[mask])
    y = np.log10(sum1[mask] / 3 + sum2[mask] / 3)

    ax2.scatter(x, y)
    ax2.set_xlabel('Sum of means (Real Data)')
    ax2.set_ylabel('Sum of means (Simulated Data)')
    slope, intercept = np.polyfit(x, y, 1)
    ax2.plot(x, x * slope + intercept, color='r', lw=3)

    return fig


def setQtrue(df, qTrue, half, nGenes):
    """ For some of the simulations I want to set 1/2 of the genes to qTrue and the other half to 1-qTrue. """
    if half:
        df.ix[:nGenes / 2 + 1, 'qTrue'] = qTrue
        df.ix[nGenes / 2 + 1:, 'qTrue'] = 1 - qTrue
    else:
        df['qTrue'] = qTrue


def main(fname, nGenes, sampleIDs, dna, bias, ai, qTrue, half, oname, figname):
    """ Main Script """
    # Open pdf
    pp = PdfPages(figname)

    # Import real data
    dat = pd.read_csv(fname)

    # Pull random subset of nGenes
    mydat = dat.iloc[np.random.choice(dat.index, nGenes, replace=False), :]
    mydat.reset_index(drop=True, inplace=True)
    del dat

    # Initialize output dataframe
    out = pd.DataFrame(columns=['fusion_id', 'qTrue', 'trueMean1', 'trueMean2'] + sampleIDs[0] + sampleIDs[1])
    out['fusion_id'] = mydat['fusion_id']
    out.reset_index(drop=True, inplace=True)

    # Calculate average across biorep for each genotype
    mean1 = mydat[sampleIDs[0]].mean(axis=1)
    mean2 = mydat[sampleIDs[1]].mean(axis=1)
    sumMean = mean1 + mean2

    # Plot original data
    mask = (mean1 > 0) & (mean2 > 0)
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    fig.suptitle('Real Data (RNA)')
    ax.plot(mean1[mask] / mean2[mask])
    ax.set_ylabel('Line / Tester before bias correction')
    ax.set_xlabel('Exonic Regions')
    pp.savefig(fig)

    # Run simulation

    if bias and ai:
        # YesBias YesAI
        R = 1.5
        setQtrue(out, qTrue, half, nGenes)

    elif bias:
        # YesBias NoAI
        R = 1
        setQtrue(out, qTrue, half, nGenes)

    else:
        # NoBias NoAI
        R = 1
        setQtrue(out, qTrue, half, nGenes=nGenes)

        ## When NoBias NoAI also generate DNA counts
        figDNA = runDNASim(out, mydat, dna, sampleIDs)
        pp.savefig(figDNA)

    figRNA = runRNASim(out, sumMean, sampleIDs, R)
    pp.savefig(figRNA)

    # Close PDF
    pp.close()

    # Write output table
    out.sort(columns='fusion_id', inplace=True)
    out.to_csv(oname, index=False)


if __name__ == '__main__':
    opts = getOptions()

    # RNA column names
    sample1 = ['line_RNA_Rep1', 'line_RNA_Rep2', 'line_RNA_Rep3']
    sample2 = ['tester_RNA_Rep1', 'tester_RNA_Rep2', 'tester_RNA_Rep3']
    sampleIDs = (sample1, sample2)

    # DNA column names
    dna1 = ['LINE_DNA_TOTAL']
    dna2 = ['TESTER_DNA_TOTAL']
    dna = (dna1, dna2)

    args = [opts.fname, opts.nGenes, sampleIDs, dna, opts.bias, opts.ai, opts.qTrue, opts.half, opts.oname, opts.figname]

    main(*args)
