#!/usr/bin/env python

import argparse
import pandas as pd

# Import custon functions for FSM reduction
import FSM_consolidation_functions as FSM

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get transcript-level fusion (exonic region) variables")

    # Input data
    parser.add_argument("-f", "--fusions", dest="inFus", required=True, help="CSV file of fusion (exonic region) annotations including transcript_id column (separated by '|' if there are multiple (*_fusion_annotations.csv)")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file of transcript-level fusion (exonic regions) variables")

    args = parser.parse_args()
    return args

def main():

    # Get input file of fusions
    fusionDF = pd.read_csv(args.inFus, low_memory=False) 

    # Add exon region length
    fusionDF['exonRegion_length'] = fusionDF['fusion_stop'] - fusionDF['fusion_start']

    # Count exon regions
    print("{} unique exon regions (fusions) present in annotation".format(fusionDF['fusion_id'].nunique()))
    
    # Split transcript_id in by pipes and keep all other values the same
    # Some fusions are found in multiple transcripts (common or constitutive) and
    #     have a piped list of transcript_id values in the transcript_id column
    splitFusionDF = FSM.split_transcript_id(df=fusionDF,sort_list=['transcript_id','fusion_start'])
    del(fusionDF)
    
    # Make transcript-level exon region variables by collapsing by transcript_id
    collapseFusionDF = FSM.collapse_fragments(df=splitFusionDF,sort_list=['transcript_id','fusion_start'],
                                              level="transcript", featureType="fusion",
                                              featureName="exonRegion").drop(columns=['gene_id'])
    del(splitFusionDF)
    
    # Output transcript-level file
    collapseFusionDF.to_csv(args.outFile, index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

