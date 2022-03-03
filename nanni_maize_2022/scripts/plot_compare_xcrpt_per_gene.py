#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Plot transcripts per gene between two datasets (reference vs. data, data1 vs. data2, etc.)")

    # Input data
    parser.add_argument("-1", "--input-file1", dest="inF1", required=True, help="CSV file for data 1 of event_id to transcript_id to gene_id (*_event2transcript2gene_index.csv)")
    parser.add_argument("-2", "--input-file2", dest="inF2", required=True, help="CSV file for data 2 of event_id to transcript_id to gene_id (*_event2transcript2gene_index.csv)")
    parser.add_argument("-n1", "--name1", dest="inName1", required=True, help="Name for data 1 to be used in title/axis")
    parser.add_argument("-n2", "--name2", dest="inName2", required=True, help="Name for data 2 to be used in title/axis")
    parser.add_argument("--hist-1-2", dest="hist12", action="store_true", required=False, help="Make histogram of transcript per gene differences between data 1 and data 2 (data 1 minus data 2)")
    parser.add_argument("--hist-2-1", dest="hist21", action="store_true", required=False, help="Make histogram of transcript per gene differences between data 2 and data 1 (data 2 minus data 1)")
    parser.add_argument("--scatter", dest="scatter", action="store_true", required=False, help="Make scatter plot of transcript per gene values of data 1 and data 2")

    # Output data
    parser.add_argument("-d", "--out-directory", dest="outDir", required=True, help="Directory for output files")
    parser.add_argument("-p", "--out-prefix", dest="outPrefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    return args

def main():
    # Get names of data and input files of event_id to transcript_id to gene_id
    name1 = args.inName1
    name2 = args.inName2
    annotDF1 = pd.read_csv(args.inF1, low_memory=False)
    annotDF2 = pd.read_csv(args.inF2, low_memory=False)

    # Extract all unique pairs of individual transcript_id to gene_id
    #     and count transcripts per gene in each dataset
    xcrptDF1 = annotDF1[(~annotDF1['transcript_id'].str.contains("|",regex=False))&
            (~annotDF1['gene_id'].str.contains("|",regex=False))&
            (annotDF1['transcript_id']!="Unannotated")][['gene_id','transcript_id']].drop_duplicates()
    geneDF1 = xcrptDF1.groupby('gene_id')['transcript_id'].count().reset_index()
    xcrptDF2 = annotDF2[(~annotDF2['transcript_id'].str.contains("|",regex=False))&
            (~annotDF2['gene_id'].str.contains("|",regex=False))&
            (annotDF2['transcript_id']!="Unannotated")][['gene_id','transcript_id']].drop_duplicates()
    geneDF2 = xcrptDF2.groupby('gene_id')['transcript_id'].count().reset_index()
    
    # Merge datasets and count
    mergeDF = pd.merge(geneDF1,geneDF2,how='outer',on='gene_id',suffixes=('_1','_2'),indicator='merge_check')
    mergeDF['transcript_id_2'] = np.where(mergeDF['merge_check']=='left_only',0,mergeDF['transcript_id_2'])
    countDF = mergeDF['merge_check'].value_counts()
    print("{} genes in {} only\n{} genes in {} only\n{} genes in both".format(
            countDF['left_only'],name1,countDF['right_only'],name2,countDF['both']))
    
    # Plot transcripts per gene values of data1 vs. data2 (genes only in one dataset will be 0 in the other)
    
    # Histogram of difference between 1 and 2
    if args.hist12:
        bins = list(range(0,1,1)) + list(range(1,22,2)) + list(range(40,120,20)) + [max(120,(mergeDF['transcript_id_1']-mergeDF['transcript_id_2']).max())]
        hist, binEdges = np.histogram((mergeDF['transcript_id_1']-mergeDF['transcript_id_2']),bins)
        fig,ax = plt.subplots()
        # Plot  histogram against ints on the x axis
        ax.bar(range(len(hist)),hist,width=1)
        # Set the ticks to the middle of the bars
        ax.set_xticks([i-0.5 for i,j in enumerate(hist)])
        # Set the xticklabels to a string for bin range
        ax.set_xticklabels(['{} - {}'.format(bins[i],bins[i+1]) for i,j in enumerate(hist)], rotation=45)
        plt.title("Difference in transcripts per gene in\n{} vs. {}".format(name1,name2))
        plt.savefig("{}/{}_hist_1-2.png".format(args.outDir,args.outPrefix),dpi=600,format="png")
        fig.clear()
    if args.hist21:
        bins = list(range(0,1,1)) + list(range(1,22,2)) + list(range(40,120,20)) + [max(120,(mergeDF['transcript_id_2']-mergeDF['transcript_id_1']).max())]
        hist, binEdges = np.histogram((mergeDF['transcript_id_2']-mergeDF['transcript_id_1']),bins)
        fig,ax = plt.subplots()
        # Plot  histogram against ints on the x axis
        ax.bar(range(len(hist)),hist,width=1)
        # Set the ticks to the middle of the bars
        ax.set_xticks([i-0.5 for i,j in enumerate(hist)])
        # Set the xticklabels to a string for bin range
        ax.set_xticklabels(['{} - {}'.format(bins[i],bins[i+1]) for i,j in enumerate(hist)], rotation=45)
        plt.title("Difference in transcripts per gene in\n{} vs. {}".format(name2,name1))
        plt.savefig("{}/{}_hist_2-1.png".format(args.outDir,args.outPrefix),dpi=600,format="png")
        fig.clear()
    if args.scatter:
        plt.scatter(x=mergeDF['transcript_id_1'],y=mergeDF['transcript_id_2'],s=10)
        plt.xlabel(name1)
        plt.ylabel(name2)
        plt.axis('image')
        plt.title("Transcripts Per Gene")
        plt.savefig("{}/{}_scatter.png".format(args.outDir,args.outPrefix),dpi=600,format="png")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
