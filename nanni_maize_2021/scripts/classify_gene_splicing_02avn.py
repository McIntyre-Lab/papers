#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Classify types of splicing in genes (alt. exon, alt. donor/acceptr, intron retention)")

    # Input data
    parser.add_argument("-d", "--distance", dest="inDist", required=True, help="Input CSV of pairwise transcript distances (*_pairwise_transcript_distance.csv)")
    parser.add_argument("-f", "--fragment", dest="inFrag", required=True, help="Input CSV of intron retention flagged fragment annotation file (*_exon_fragment_annotations_flag_IR.csv)")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file")

    args = parser.parse_args()
    return args

def main():
    # Get input pairwise distance file and check if it is empty
    distDF = pd.read_csv(args.inDist, low_memory=False)
    if len(distDF) == 0:
        print("ERROR: Empty distance file provided ... No plots generated")
        sys.exit()
        
    # Get number of transcripts per gene from distance file
    # Number of pairwise comparisons p is p = k(k-1)/2 where k is the number of transcripts
    #     so k = (1 + sqrt(1+8p))/2
    distDF['numPair'] = distDF.groupby('gene_id')['transcript_1'].transform('count')
    distDF['numXcrpt'] = (1+np.sqrt(1+(8*distDF['numPair'])))/2
    del(distDF['numPair'])

    # Get input fragment annotation file flagged with involvement in IR
    fragDF = pd.read_csv(args.inFrag, low_memory=False)
    
    # Flag NIC/NNC transcripts ("PB" in transcript_id)
    fragDF['flag_has_NIC_NNC'] = np.where(fragDF['transcript_id'].str.contains("PB"),1,0)
    
    # Flag if NIC/NNC has putative IR fusion
    fragDF['flag_NIC_NNC_IR_fusion'] = np.where((fragDF['flag_has_NIC_NNC']==1)&
          (fragDF['flag_intron_retention']==1),1,0)
    
    # Get unique gene_id values with IR flags
    geneFragDF = fragDF[fragDF['flag_multigene']==0].groupby('gene_id')[['flag_intron_retention',
                               'flag_has_NIC_NNC','flag_NIC_NNC_IR_fusion']].max().reset_index()
    
    # Flag alternative exons (prop_ER_diff > 0)
    distDF['flag_alt_exon'] = np.where(distDF['prop_ER_diff']>0,1,0)
    
    # Flag alternative donor/acceptors (num_nt_diff_in_shared_ER > 0)
    distDF['flag_alt_donor_acceptor'] = np.where(distDF['num_nt_diff_in_shared_ER']>0,1,0)
    
    # Get unique gene_id values with alt exon/donor/acceptor flags
    geneDistDF = distDF.groupby('gene_id')[['flag_alt_exon','flag_alt_donor_acceptor',
                               'numXcrpt','num_nt_diff_in_shared_ER']].max().reset_index().rename(columns={'num_nt_diff_in_shared_ER':'max_num_nt_diff_in_shared_ER'})
    
    # Merge gene dataframes on gene_id and only keep those in distance file (have >=2 transcripts)
    mergeDF = pd.merge(geneFragDF,geneDistDF,how='outer',on='gene_id',indicator='merge_check')
    mergeDF = mergeDF[mergeDF['merge_check']=="both"]
    del(mergeDF['merge_check'])
    
    # Check that no genes have none (all should have at least one type)
    if len(mergeDF[(mergeDF['flag_intron_retention']==0)&(mergeDF['flag_alt_exon']==0)&(mergeDF['flag_alt_donor_acceptor']==0)]) > 0:
        print("ERROR: At least one gene with no splicing classifications found")
        sys.exit()
    
    # Output csv file of flags
    mergeDF.to_csv(args.outFile,index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
