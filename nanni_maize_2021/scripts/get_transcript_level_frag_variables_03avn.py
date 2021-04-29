#!/usr/bin/env python

import argparse
import pandas as pd

# Import custon functions for FSM reduction
import FSM_consolidation_functions as FSM

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get transcript-level fragment variables")

    # Input data
    parser.add_argument("-f", "--fragments", dest="inFrag", required=True, help="CSV file of fragments annotations including transcript_id column (separated by '|' if there are multiple (*_exon_fragment_annotations.csv)")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file of transcript-level fragment variables")

    args = parser.parse_args()
    return args

def main():
    # Get input file of fragments
    fragDF = pd.read_csv(args.inFrag, low_memory=False)
    
    # Add fragment length
    fragDF['fragment_length'] = fragDF['fragment_stop'] - fragDF['fragment_start']
    
    # Count fragments
    print("{} unique exon fragments present in annotation".format(fragDF['fragment_id'].nunique()))
    
    # Split transcript_id in by pipes and keep all other values the same
    # Some fragments are found in multiple transcripts (common or constitutive) and
    #     have a piped list of transcript_id values in the transcript_id column
    splitFragDF = FSM.split_transcript_id(df=fragDF,sort_list=['transcript_id','fragment_start'])
    del(fragDF)
    
    # Make transcript-level fragment variables by collapsing by transcript_id
    collapseFragDF = FSM.collapse_fragments(df=splitFragDF,sort_list=['transcript_id','fragment_start'],
                                            level="transcript").drop(columns=['gene_id'])
    del(splitFragDF)
    
    # Output transcript-level file
    collapseFragDF.to_csv(args.outFile, index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

