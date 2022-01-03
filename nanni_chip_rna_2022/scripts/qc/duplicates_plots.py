#!/usr/bin/env python

##        DESCRIPTION: This program generates plots for assessing the proportion of unique reads in a set of
##        FASTQ files. It takes a CSV of percentage of unique reads per sample and will output two plots: a bar chart 
##        showing the percentage of unique reads by sample and a box plot that looks at the distribution of 
##        unique read content across the set of FASTQ files.
##        
##        Author: Jeremy R. B. Newman


## Import packages
import argparse
import numpy as np
import scipy
import pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import CSV and give output names")
    parser.add_argument("--input", dest="csvInput", action='store', required=True, help="Input CSV")
    parser.add_argument("--obar", dest="outBar", action='store', required=True, help="Output filename of bar chart")
    parser.add_argument("--ohist", dest="outHisto", action='store', required=True, help="Output filename of histogram")
    parser.add_argument("--name", dest="setName", action='store', required=True, help='Name for the set of FASTQ files (e.g. "Plate 2")')

    args = parser.parse_args()
    return args

def main():
    # import data
    duplicateSummary=pandas.read_csv(args.csvInput, sep=',', header=None, skiprows=1)
    
    # Plot sequence uniqueness by sample
    xLabels=duplicateSummary[0]
    percUniqueSeq=duplicateSummary[3]
    N=len(percUniqueSeq)
    
    ind = np.arange(N)
    width = 0.1
   
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    
    width = .25
    ind = np.arange(N) 
    plt.bar(ind, percUniqueSeq, width=width, color='k')
    plt.xticks(ind + width / 2, xLabels, rotation='vertical', size=6)
    barTitle='Read uniqueness by sample from ' + str(args.setName)
    ax.set_title(barTitle)
    
    ax.set_xlabel('SampleID')
    ax.set_ylabel('Unique reads (%)')
    plt.ylim(40,100)
    
    fig.savefig(args.outBar)

    # Plot the distribution of sequence uniqueness across all samples
    percUniqueSeqDist=duplicateSummary[3].tolist()
    
    # Create a figure instance
    fig1 = plt.figure(figsize=(9, 6))
    
    # Create an axes instance
    ax1 = fig1.add_subplot(111)

    ## add patch_artist=True option to ax.boxplot() 
    ## to get fill color
    bp = ax1.boxplot(percUniqueSeqDist, patch_artist=True)
    
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#000000', linewidth=2)
        # change fill color
        box.set( facecolor = '#fc5050' )
    
    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#000000', linewidth=2, linestyle='solid')
    
    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#000000', linewidth=2)
    
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#000000', linewidth=2)
    
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#0000ff', alpha=0.75)
        
    histoTitle='Distribution of read uniqueness across ' + str(args.setName)
    plt.title(histoTitle)
    
    plt.ylabel('Unique reads (%)')
    plt.ylim(40,100)
    plt.xlabel('')
    fig1.savefig(args.outHisto)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

