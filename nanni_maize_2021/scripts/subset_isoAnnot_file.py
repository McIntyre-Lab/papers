#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Subset previous isoAnnot file for maize PacBio transcriptome")

    # Input data
    parser.add_argument("-i", "--isoAnnot", dest="inAnnot", required=True, help="Input GFF-like file of functional annotations for maize PacBio transcriptome (Has all PB transcript id values)")
    parser.add_argument("-c", "--class", dest="inClass", required=True, help="Input SQANTI QC classification file for maize PacBio transcriptome (Has all PB transcript id values that are in isoAnnot input file)")
    parser.add_argument("-l", "--list", dest="inList", required=True, help="Input list of transcript_id values for consolidated maize PacBio transcriptome (Has some reference transcript_id values for FSM/ISM isoforms that associated with them)")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output GFF-like file of functional annotations for consolidated maize PacBio transcriptome (with PB transcript id values changed to associated reference when appropriate)")

    args = parser.parse_args()
    return args

def main():
    # Get input isoAnnot file
    annotDF = pd.read_csv(args.inAnnot, sep="\t",
                          names=['isoform','source','feature','start','end','score',
                                 'strand','frame','attribute'],low_memory=False)
    
    # Get input classification file
    classDF = pd.read_csv(args.inClass, sep="\t")[['isoform','structural_category',
                         'associated_gene','associated_transcript']]
    
    # Get list of transcripts to include in subset isoAnnot file
    listDF = pd.read_csv(args.inList,names=['transcript_id'])
    
    # Merge isoAnnot and classification by isoform
    mergeClassDF = pd.merge(annotDF,classDF,how='outer',on='isoform',indicator='merge_check')
    
    # Get new isoform ID by substituting FSM/ISM PB id with associated transcript id
    mergeClassDF['new_isoform'] = np.where(mergeClassDF['associated_transcript']=="novel",
                mergeClassDF['isoform'],mergeClassDF['associated_transcript'])
    
    # Merge list of transcript to include in subset with isoAnnot by new_isoform and transcript_id
    mergeAllDF = pd.merge(mergeClassDF,listDF,how='outer',left_on='new_isoform',
                          right_on='transcript_id',indicator='merge_check2')
    # mergeAllDF['merge_check2'].value_counts()
    # Checked counts and there are no 'right_only' meaning all transcripts in list were found
    
    # Output subset isoAnnot with new isoform ID
    mergeAllDF[mergeAllDF['merge_check2']=='both'][['new_isoform','source',
              'feature','start','end','score','strand','frame','attribute']].to_csv(
              args.outFile,sep="\t",index=False, header=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()