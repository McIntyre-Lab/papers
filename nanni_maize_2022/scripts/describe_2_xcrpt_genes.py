#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Describe the splicing of genes with 2 transcripts")

    # Input data
    parser.add_argument("-f", "--flags", dest="inFlag", required=True, help="Input CSV of splicing type flags")

    # Output data
#    parser.add_argument("-", "--", dest="", required=True, help="")

    args = parser.parse_args()
    return args

def main():
    # Get input files
    flagDF = pd.read_csv(args.inFlag, low_memory=False)
    flagDF[[c for c in flagDF.columns if "flag" in c]] = flagDF[[c for c in flagDF.columns if "flag" in c]].astype(int)
    
    # Output groups of splicing flags from genes with 2 transcripts
    print("{} out of {} genes with two transcripts have at least one NIC/NNC".format(
            len(flagDF[(flagDF['flag_has_NIC_NNC']==1)&(flagDF['numXcrpt']==2)]),
            len(flagDF[flagDF['numXcrpt']==2])))
    twoGene = flagDF[flagDF['numXcrpt']==2].groupby(['flag_intron_retention',
                    'flag_alt_exon','flag_alt_donor_acceptor','flag_has_NIC_NNC',
                    'flag_NIC_NNC_IR_fusion'])['gene_id'].count()
    twoGene = twoGene.append(pd.Series({('total','','','',''):len(flagDF[flagDF['numXcrpt']==2])}))
    print("Groups of splicing flags from genes with 2 transcripts:\n\n{}".format(
            twoGene.reset_index().to_string(index=False)))

    # Output distribution of nt difference in genes with only alt donor/acceptor
    print("Distribution of nt differences (in shared ER) in all 2 transcript genes with only alt. donor/acceptor:\n{}".format(
            flagDF[(flagDF['numXcrpt']==2)&(flagDF['flag_intron_retention']==0)&
                   (flagDF['flag_alt_exon']==0)&(flagDF['flag_alt_donor_acceptor']==1)]['max_num_nt_diff_in_shared_ER'].describe().to_string()))
    print("Distribution of nt differences (in shared ER) in all 2 transcript genes with only alt. donor/acceptor (NIC/NNC only):\n{}".format(
            flagDF[(flagDF['numXcrpt']==2)&(flagDF['flag_intron_retention']==0)&
                   (flagDF['flag_alt_exon']==0)&(flagDF['flag_alt_donor_acceptor']==1)&
                   (flagDF['flag_has_NIC_NNC']==1)]['max_num_nt_diff_in_shared_ER'].describe().to_string()))
    
    # Count number of genes with different splicing types in groups of transcripts per gene
#    groupDF = flagDF.groupby('numXcrpt').agg({'gene_id':'count',
#                            'flag_intron_retention':'sum',
#                            'flag_alt_exon':'sum',
#                            'flag_alt_donor_acceptor':'sum'})

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
