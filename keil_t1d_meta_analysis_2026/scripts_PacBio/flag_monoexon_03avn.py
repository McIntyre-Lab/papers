#!/usr/bin/env python
 
import argparse
import pandas as pd
import numpy as np
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Flag monoexon PacBio isoforms and reference transcripts in and output counts")

    # Input arguments
    parser.add_argument("-r", "--reference", dest="inRef", required=True, help="Reference GTF file to get reference exon counts")
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Input GTF file to make monoexon flags for")

    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output flag files")

    args = parser.parse_args()
    return args


def main():
    # Check that input files exits
    if not os.path.isfile(args.inFile) or not os.path.isfile(args.inRef):
        raise FileNotFoundError

    # Get exon counts from reference GTF
    refgtf = pd.read_csv(args.inRef,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    refexon = refgtf[refgtf["feature"]=="exon"].copy()
#    print("Total lines in reference GTF= "+str(len(gtf)))

    # Get gene_id and transcript_id values from attributes column
    for i in refexon.index:
        raw_attrs = refexon.at[i, 'attribute']
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = np.nan, np.nan
        for item in g_t_attrs:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
            elif 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
        if transcript_id == np.nan and refexon.at[i, 'feature'] != "gene" :
            print("WARNING: transcript_id not found in {}".format(refexon[i]))
        refexon.at[i, 'gene_id'] = str(gene_id)
        refexon.at[i, 'transcript_id'] = str(transcript_id)
#    print("Total transcripts in original GTF= "+str(refexon["transcript_id"].nunique()))
#    print("Total genes in original GTF= "+str(refexon["gene_id"].nunique()))

    # Get min exons per transcript for each gene to flag genes with at least one monoexon transcript
    refxcrpt = refexon.groupby(["gene_id", "transcript_id"])["feature"].count().reset_index().rename(columns={"feature": "num_exon"})
    refgene = refxcrpt.groupby("gene_id")["num_exon"].min().reset_index().rename(columns={"num_exon": "min_exon"})
    refgene["flag_ref_monoexon"] = np.where(
        refgene["min_exon"] == 1,
        1,
        0
    )

    # Get input GTF file
    ingtf = pd.read_csv(args.inFile,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    inexon = ingtf[ingtf["feature"]=="exon"].copy()
#    print("Total lines in reference GTF= "+str(len(gtf)))

    # Get gene_id and transcript_id values from attributes column
    for i in inexon.index:
        raw_attrs = inexon.at[i, 'attribute']
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = np.nan, np.nan
        for item in g_t_attrs:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
            elif 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
        if transcript_id == np.nan and inexon.at[i, 'feature'] != "gene" :
            print("WARNING: transcript_id not found in {}".format(inexon[i]))
        inexon.at[i, 'gene_id'] = str(gene_id)
        inexon.at[i, 'transcript_id'] = str(transcript_id)
#    print("Total transcripts in original GTF= "+str(inexon["transcript_id"].nunique()))
#    print("Total genes in original GTF= "+str(inexon["gene_id"].nunique()))

    # Get exons per transcript
    inxcrpt = inexon.groupby(["gene_id", "transcript_id"])["feature"].count().reset_index().rename(columns={"feature": "num_exon"})
    inxcrpt["flag_PB_monoexon"] = np.where(
        inxcrpt["num_exon"] == 1,
        1,
        0
    )

    # Merge input with reference genes
    mergeRefIn = pd.merge(
        refgene,
        inxcrpt,
        how="outer",
        on="gene_id",
        indicator="merge_check",
        validate="1:m"
    )
    # Flag unspliced fragments, where ref is multi-exon and input transcript in monoexon
    outDF = mergeRefIn[mergeRefIn["merge_check"]!="left_only"].drop(columns=["merge_check"])
    outDF["flag_unspliced_fragment"] = np.where(
        (outDF["flag_PB_monoexon"] == 1)
        & (outDF["flag_ref_monoexon"] == 0),
        1,
        0
    )
    outDF.to_csv(args.outFile, index=False)
    
     
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

