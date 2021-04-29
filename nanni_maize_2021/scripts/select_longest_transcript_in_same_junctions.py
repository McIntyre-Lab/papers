#!/usr/bin/env python

import argparse
import pandas as pd

def restricted_float(val):
    try:
        val = float(val)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (val,))

    if val <= 0.0 or val >= 1.0:
        raise argparse.ArgumentTypeError("%r not in range (0.0, 1.0)" % (val,))
    return val

def restricted_int(val):
    try:
        val = int(val)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not an integer literal" % (val,))
        
    if val < 0:
        raise argparse.ArgumentTypeError("%r not positive integer" % (val,))
    return val

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Select longest transcript within groups of transcripts with identical junctions, < N nt different, and < P% of total nt different")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Input CSV of pairwise transcript distances (*_pairwise_transcript_distance.csv)")
    parser.add_argument("-N", "--max-nt-diff", dest="inNT", required=False, type=restricted_int, help="Maximum number of nucleotides different (non-inclusive, < N) between transcripts")
    parser.add_argument("-P", "--max-percent-nt-diff", dest="inPercent", type=restricted_float, required=False, help="Maximum percent of total nucleotides different (non-inclusive, < N%) between transcripts, value must be in range (0.0, 1.0)")

    # Output data
    parser.add_argument("-d", "--output-directory", dest="outDir", required=True, help="Output directory")
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    return args

def main():
    # Get input pairwise distance file
    distDF = pd.read_csv(args.inFile, low_memory=False)
    
    # Get threshold values
    maxNT = args.inNT
    if args.inPercent is not None:
        maxPercDiff = args.inPercent
        maxPercSim = 1 - maxPercDiff
    else:
        maxPercDiff = None
        maxPercSim = None
        
    # Calculate number of nucleotides different for each pair
    distDF['num_nt_diff'] = distDF['num_nt_T1_only'] + distDF['num_nt_T2_only']
    distDF['num_nt_diff_in_shared_ER'] = distDF['num_nt_T1_only_in_shared_ER'] + distDF['num_nt_T2_only_in_shared_ER']
    
    # Select transcript pairs with:
    # 1) Identical junctions (prop_junction_similar == 1)
    # 2) if set, < N nt different (num_nt_diff < N)
    # 3) if set, < P% of total nt different (prop_nt_similar >= (1-P))
    subsetDF = distDF[distDF['prop_junction_similar']==1]
    if maxNT is not None:
        subsetDF = subsetDF[subsetDF['num_nt_diff']<maxNT]
    if maxPercSim is not None:
        subsetDF = subsetDF[subsetDF['prop_nt_similar']>=maxPercSim]

    # Get individual transcripts from the pairs
    t1DF = subsetDF[['gene_id','transcript_1','junction_shared','num_nt_shared','num_nt_T1_only']].copy().rename(columns={'transcript_1':'transcript_id'})
    t2DF = subsetDF[['gene_id','transcript_2','junction_shared','num_nt_shared','num_nt_T2_only']].copy().rename(columns={'transcript_2':'transcript_id'})
    t1DF['length'] = t1DF['num_nt_shared'] + t1DF['num_nt_T1_only']
    t2DF['length'] = t2DF['num_nt_shared'] + t2DF['num_nt_T2_only']
    del(t1DF['num_nt_T1_only'])
    del(t1DF['num_nt_shared'])
    del(t2DF['num_nt_T2_only'])
    del(t2DF['num_nt_shared'])

    # Concatenate and get unique transcripts
    concatDF = pd.concat([t1DF,t2DF],ignore_index=True).drop_duplicates()
    
    # Add list of transcripts with shared ER sets
    concatDF = concatDF.sort_values(['gene_id','junction_shared','length'],ascending=False)
    concatDF['transcript_set'] = concatDF.groupby(['gene_id','junction_shared'])['transcript_id'].transform(lambda x: "|".join(x))
    
    # Groupby gene_id and set of ER then select the longest representative
    longestDF = concatDF.loc[concatDF.groupby(['gene_id','junction_shared'])['length'].idxmax()]
    
    # Output counts
    print("\n{} input transcripts in pairwise comparisons in {} genes (genes must have at least 2 transcripts)".format(distDF['transcript_1'].append(distDF['transcript_2']).nunique(),distDF['gene_id'].nunique()))
    print("{} transcripts in {} genes that share identical junction sets with at least one other transcripts and fall within thresholds (if set)".format(concatDF['transcript_id'].nunique(),concatDF['gene_id'].nunique()))
    print("{} longest representative transcripts selected from the {} transcripts".format(len(longestDF),concatDF['transcript_id'].nunique()))
    print("{} of the {} multi-transcript genes now have 1 representative transcript".format(longestDF[~longestDF['gene_id'].isin(distDF[distDF['prop_junction_similar']!=1]['gene_id'])]['gene_id'].nunique(),distDF['gene_id'].nunique()))
    
    # Output longest transcript_id to list of transcripts
    suffix = ""
    if maxNT is not None:
        suffix = "_{}LT{}nt".format(suffix,maxNT)
    if maxPercDiff is not None:
        suffix = "_{}LT{}percNT".format(suffix,int(maxPercDiff*100))
    longestDF[['gene_id','transcript_id','transcript_set']].to_csv(
            "{}/{}_transcriptID_2_transcriptSet_sharedER{}.csv".format(args.outDir,args.outPrefix,suffix),index=False)
    
    # Output list of transcript_id values dropped (were not the longest)
    notRetained = concatDF[~concatDF['transcript_id'].isin(longestDF['transcript_id'])]['transcript_id']
    notRetained.to_csv("{}/{}_transcript_list_not_retained.csv".format(args.outDir,args.outPrefix),index=False,header=False)
    
    # Output reduced set of transcript distance
    reducedDF = distDF[(~distDF['transcript_1'].isin(notRetained))&
                       (~distDF['transcript_2'].isin(notRetained))]
    reducedDF.to_csv("{}/{}_pairwise_transcript_distance_reduced_sharedJunction.csv".format(args.outDir,args.outPrefix),index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
