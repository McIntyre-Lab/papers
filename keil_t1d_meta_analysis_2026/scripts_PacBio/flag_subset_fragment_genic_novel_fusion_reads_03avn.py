#!/usr/bin/env python
 
import argparse
import pandas as pd
import numpy as np
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
                "Flag reads in putative novel/fusion loci or are "
                "genic/genic_intron, flag reads that pass selection, "
                "and output counts and list of passing transcripts."
        )
    )

    # Input arguments
    parser.add_argument(
        "-c",
        "--classification",
        dest="inFile",
        required=True,
        help="Monoexon flagged classification file from SQANTI QC"
    )
    parser.add_argument(
        "-n",
        "--keep-novel",
        dest="keep_novel",
        action="store_true",
        help=(
            "Keep reads from novel loci in selection. "
            "Default: Remove novel loci associated reads"
        )
    )
    parser.add_argument(
        "-f",
        "--keep-fusion",
        dest="keep_fusion",
        action="store_true",
        help=(
            "Keep reads from fusion loci in selection (NOTE: fusion loci are "
            "identified by the presence of an underscore in the associated "
            "gene, if reference gene_id values also contain underscores then "
            "this identification of fusion loci will be inaccurate). "
            "Default: Remove fusion loci associated reads"
        )
    )
    parser.add_argument(
        "-g",
        "--keep-genic",
        dest="keep_genic",
        action="store_true",
        help=(
            "Keep reads classified as genic/genic_intron in selection. "
            "Default: Remove genic/genic_intron reads"
        )
    )
    parser.add_argument(
        "-u",
        "--keep-unspliced",
        dest="keep_unspliced",
        action="store_true",
        help=(
            "Keep reads classified as unspliced fragments (or monoexon reads "
            "associated with a multi-exon reference) in selection. "
            "Default: Remove genic/genic_intron reads"
        )
    )

    # Output arguments
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output file for frequencies"
    )
    parser.add_argument(
        "-p",
        "--output_reads_passed",
        dest="outPassed",
        required=True,
        help="Output file for reads passing selection parameters"
    )

    args = parser.parse_args()
    return args

def rchop(string, ending):
    if string.endswith(ending):
        return string[:-len(ending)]
    return string

def main():
    # Check that input files exits
    if not os.path.isfile(args.inFile) :
        raise FileNotFoundError

    # Get input monoexon flagged classification file
    classDF = pd.read_csv(args.inFile, sep="\t", low_memory=False)
    
    # Flag unspliced fragments
    # Where the read is monoexon and the associated reference is multiexon
    classDF['flag_unspliced_fragment'] = np.where(
            (classDF['flag_PB_monoexon']==0)|
                    (classDF['flag_ref_monoexon']==1),
            0,
            1
    )
    
    # Flag novel genes (genes that were annotated by SQANTI3 QC as "novelGene")
    #     and fusion genes (genes that were annotated by SQANTI3 QC with [gene1]_[gene2])
    # !!!! Fusion gene ID's are found assuming that reference gene_id values do not contain underscores
    #         (this is true for Maize B73 and C. elegans WBcel235 reference annotations)
    # Use novel_in_catalog to determine if reference gene_id values contain "_"
    if len(classDF[classDF['structural_category'].isin([
            'novel_in_catalog', 'novel_not_in_catalog'])]) > 0:
        # Check if underscores in reference gene_id values
        if len(classDF[(classDF['structural_category'].isin([
                'novel_in_catalog', 'novel_not_in_catalog']))&
                (classDF['associated_gene'].str.contains("_"))&
                (~classDF['associated_gene'].str.contains("novel"))]) > 0:
            print("!!!WARNING: Reference gene_id values contain underscores. "
                  "Identification of fusion loci in non-fusion reads "
                  "(FSM or ISM) is not possible and will be only flagged for "
                  "reads identified as fusion in SQANTI3 structural_category."
            )
            classDF['flag_associated_gene_fusion'] = np.where(
                    classDF['structural_category']=="fusion",
                    1,
                    0
            )
        else:
            classDF['flag_associated_gene_fusion'] = np.where(
                    (~classDF['associated_gene'].str.contains("novelGene"))&
                        (classDF['associated_gene'].str.contains("_")),
                    1,
                    0
            )
    classDF['flag_associated_gene_novel'] = np.where(classDF['associated_gene'].str.contains("novelGene"),1,0)
    classDF['flag_associated_gene_novel_or_fusion'] = np.where((classDF['flag_associated_gene_novel']==1)|(classDF['flag_associated_gene_fusion']==1),1,0)

    # Flag "genic" or "genic intron" reads (these are reads that do not match any junction/donor/acceptor of the reference)
    classDF['flag_genic_read'] = np.where(
            (classDF['structural_category']=="genic")|
                (classDF['structural_category']=="genic_intron"),
                1,
                0
    )
    classDF['flag_genic_read_no_fusion'] = np.where(
            (
                (classDF['structural_category']=="genic")|
                (classDF['structural_category']=="genic_intron")
            )&
            (classDF['flag_associated_gene_fusion']==0),
                1,
                0
    )

    # Flag reads that pass baseline selection
    remove_col = []
    if not args.keep_novel:
        remove_col.append('flag_associated_gene_novel')
    if not args.keep_fusion:
        remove_col.append('flag_associated_gene_fusion')
    if not args.keep_genic:
        remove_col.append('flag_genic_read')
    if not args.keep_unspliced:
        remove_col.append('flag_unspliced_fragment')
    classDF['flag_read_pass_baseline_selection'] = np.where(
            classDF[remove_col].sum(axis=1)==0,
            1,
            0
    )

    # Get frequencies of flags
    f = open(args.outFile,'w')
    f.write("{}\n\n".format(
        classDF.groupby([
            'flag_unspliced_fragment',
            'flag_associated_gene_novel',
            'flag_associated_gene_fusion',
            'flag_genic_read',
            'flag_read_pass_baseline_selection'
        ])['isoform'].count().reset_index().rename(columns={
            'isoform':'num_reads'
    }).to_string(index=False)
    ))
    f.write("{} total reads\n\n".format(
            len(classDF)))
    f.write((
        "\n{} reads pass selection\nCounts of classification categories "
        "of reads that pass selection:\n{}\n"
        ).format(
            len(classDF[classDF['flag_read_pass_baseline_selection']==1]),
            classDF[classDF['flag_read_pass_baseline_selection']==1][
                    'structural_category'].value_counts().to_string())
    )
    f.close()
    
    # Output flag classification file
    prefix = rchop(rchop(args.inFile,"flag_mono.tsv"),"flag_mono.txt")
    classDF.to_csv(prefix+"flag_monoexon_novel_fusion_loci_genic_reads.txt",
                   sep='\t',
                   index=False)
    
    # Output reads passing baseline selection
    classDF[classDF['flag_read_pass_baseline_selection']==1]['isoform'].to_csv(
            args.outPassed, index=False, header=False)
    
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

