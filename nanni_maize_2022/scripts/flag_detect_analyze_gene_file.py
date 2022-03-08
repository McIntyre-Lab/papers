#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
#from matplotlib import pyplot as plt
#from upsetplot import UpSet

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Flag maize genotype detection (CCS gene read count > 0 and/or "
            "RNA-seq gene TPM > 0 in at least one sample) and analyzability "
            "(RNA-seq TPM>5 in at least 50% of replicates in at least one "
            "transcript of the gene)."
        )
    )

    # Input data
    parser.add_argument(
        "-f",
        "--flag",
        dest="inFlag",
        required=True,
        help="Input gene flag file that includes PAV/paralog flags from Hoopes."
    )
    parser.add_argument(
        "-m",
        "--mean-data",
        dest="inMean",
        required=True,
        help=(
            "Input CSV of means calculated from tappAS with gene 'detection' "
            "flags that now correspond to analyzable flags (TPM>5 in at least "
            "50% of replicates in at least one transcript of the gene)."
        )
    )
    parser.add_argument(
        "-c",
        "--ccs-reads",
        dest="inCCS",
        required=True,
        help=(
            "Input TSV file of CCS HTSeq counts on genes of the assembled "
            "transcriptome."
        )
    )
    parser.add_argument(
        "-s",
        "--short-reads",
        dest="inShort",
        required=True,
        help=(
            "Input TSV file of short read HTSeq counts on genes of the assembled "
            "transcriptome."
        )
    )
    parser.add_argument(
        "-g",
        "--gtf",
        dest="inGTF",
        required=True,
        help=(
            "Input GTF file of assembled transcriptome (12604 gnees)."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output CSV file for all flags and means."
    )

    args = parser.parse_args()
    return args


def get_gene_from_gtf(infile):
    # Get input GTF file
    gtf = pd.read_csv(
            infile,
            names=[
                "chr",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute"
            ],
            comment="#",
            dtype=str,
            compression="infer",
            sep="\t",
            low_memory=False
    )
    # Select only exon features
    gtf = gtf[gtf["feature"]=="exon"]
    # Get attribute values
    gtf["attribute_values"] = gtf["attribute"].str.split(" ")
    # Check that gene_id and transcript_id attributes are present for all exons
    if not gtf["attribute_values"].apply(lambda x: "transcript_id" in x).all():
        print(
            "ERROR: transcript_id not contained within all exon "
            "attributes".format(infile)
        )
    if not gtf["attribute_values"].apply(lambda x: "gene_id" in x).all():
        print(
            "ERROR: gene_id not contained within all exon "
            "attributes".format(infile)
        )
    # Extract transcript_id and gene_id from attribute column
    gtf["transcript_id"] = gtf.apply(
            lambda x: x["attribute_values"][
                    x["attribute_values"].index("transcript_id")+1
                    ].split(";")[0].split("\"")[1],
            axis=1
    )
    gtf["gene_id"] = gtf.apply(
            lambda x: x["attribute_values"][
                    x["attribute_values"].index("gene_id")+1
                    ].split(";")[0].split("\"")[1],
            axis=1
    )
    # sort by gene and transcript id
    gtf = gtf.sort_values(["gene_id", "transcript_id"])
    # select first transcript_id per gene
    geneDF = gtf.groupby("gene_id")[
        [
            "chr",
            "transcript_id"
        ]].first().reset_index().rename(columns={
            "transcript_id": "first_transcript_id"
    })
    # Return unique list of genes with corresponding chromosome and first transcript
    return geneDF


def main():
    # Get input files
    flagDF = pd.read_csv(args.inFlag, low_memory=False)
    # make_combination_flag_file/combination_flag_file_shrtRd_ccs_hoopes.csv
    meanDF = pd.read_csv(args.inMean, low_memory=False)
    # 2018/PacBio/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/tappas_output_TMM_norm_1CPMfilter/plot_maize_PB_groups_UpSet_v3/zmtr_fsm_ism_nic_nnc_consol_geneCollapse_from_meanDF.csv
    ccsDF = pd.read_csv(args.inCCS, sep="\t", low_memory=False)
    # 2018/PacBio/align_raw_reads_2_genomes/htseq_gene_counts/all_ccs_gene_htseq_count.tsv
    shortDF = pd.read_csv(args.inShort, sep="\t", low_memory=False)
    # 2018/PacBio/align_raw_reads_2_genomes/htseq_gene_counts/all_shrtRd_gene_htseq_count.tsv
    assembGeneDF = get_gene_from_gtf(args.inGTF)
    # 2018/PacBio/sqanti_classification_category_subset/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf

    # Select only PAV and paralog flags/variables from flag file
    flagDF = flagDF[[
            "geneID",
            "flag_zea_mays_paralog_hoopes",
            "orthogroup_hoopes",
            "flag_pav",
            "freq_pav_hoopes"
    ]]

    # Merge PAV/paralog flags with means
    assembFlagMerge = pd.merge(
            flagDF,
            assembGeneDF,
            how="outer",
            left_on="geneID",
            right_on="gene_id",
            indicator="merge_check",
    )
    # Check merge
    if assembFlagMerge["merge_check"].value_counts()["both"] != len(assembGeneDF):
        print("ERROR: Not all genes in GTF are found in flag file.")
    # Drop genes not in the 12604
    assembFlagDF = assembFlagMerge[
            assembFlagMerge["merge_check"]=="both"
    ].drop(columns=["merge_check", "geneID"])

    for genotype in ["B73", "C123", "Hp301", "Mo17", "NC338"]:
        # Rename flag_detect_ to flag_analyze_ in the mean file (this flag
        #    represents genes with TPM>5 in at least 50% of replicates in at
        #    least one transcript)
        meanDF = meanDF.rename(columns={
                "flag_detect_" + genotype: "flag_analyze_" + genotype
        })

        # Make RSEM (to transcripts) shrtRd detection flags for genes with mean TPM>0 in either treatment
        meanDF["flag_detect_shrtRd_rsem_" + genotype] = np.where(
            (meanDF["mean_" + genotype + "_Amb"] +
                meanDF["mean_" + genotype + "_Ele"] > 0),
            1,
            0
        )

    meanDF = meanDF.rename(columns={"sum_detect": "sum_analyze"})
    
    # Merge means with assembled transcript gene flags
    assembFlagMeanMerge = pd.merge(
            assembFlagDF,
            meanDF,
            how="outer",
            on="gene_id",
            indicator="merge_check",
            validate="1:1"
    )
    # Check merge
    if assembFlagMeanMerge["merge_check"].value_counts()["both"] != len(meanDF):
        print("ERROR: Not all genes in GTF are found in mean file.")
    # Drop merge_check columns
    assembFlagMeanMerge = assembFlagMeanMerge.drop(columns=["merge_check"])

    # Flag detection of ccs flags in htseq output if the sum of the amb and oz
    #     samples in the genotype are > 0 (at least one read)
    for genotype in ["B73", "C123", "Hp301", "Mo17", "NC338"]:
        ccsDF["flag_detect_ccs_" + genotype] = np.where(
            (ccsDF[[
                    c for c in ccsDF.columns if genotype.lower() in c
                ]].sum(axis=1) > 0),
            1,
            0
        )

    # Merge ccs detection flags with assembled transcript gene flags and means
    assembFlagMeanCCS = pd.merge(
            assembFlagMeanMerge,
            ccsDF,
            how="outer",
            on="gene_id",
            indicator="merge_check",
            validate="1:1"
    )
    # Check merge
    if assembFlagMeanCCS["merge_check"].value_counts()["both"] != len(assembFlagMeanMerge):
        print("ERROR: Not all genes in GTF are found in ccs file.")
    # Drop merge_check columns
    assembFlagMeanCCS = assembFlagMeanCCS[
            assembFlagMeanCCS["merge_check"] == "both"
        ].drop(columns=["merge_check"])

    # Flag detection of short read flags in htseq (genome alignment) output if the sum of all amb and oz
    #     replicate samples in the genotype are > 0 (at least one read)
    for genotype in ["B73", "C123", "Hp301", "Mo17", "NC338"]:
        shortDF["flag_detect_shrtRd_genomeAln_" + genotype] = np.where(
            (shortDF[[
                    c for c in shortDF.columns if genotype in c
                ]].sum(axis=1) > 0),
            1,
            0
        )

    # Merge short read detection flags with assembled transcript gene flags, means, and ccs detection flags
    assembFlagMeanDetect = pd.merge(
            assembFlagMeanCCS,
            shortDF,
            how="outer",
            on="gene_id",
            indicator="merge_check",
            validate="1:1"
    )
    # Check merge
    if assembFlagMeanDetect["merge_check"].value_counts()["both"] != len(assembFlagMeanCCS):
        print("ERROR: Not all genes in GTF are found in ccs file.")
    # Drop merge_check columns
    assembFlagMeanDetect = assembFlagMeanDetect[
            assembFlagMeanDetect["merge_check"] == "both"
        ].drop(columns=["merge_check"])

    # Add genotype detection flags using both ccs and short reads
    detectFlagCols = []
    for genotype in ["B73", "C123", "Hp301", "Mo17", "NC338"]:
        assembFlagMeanDetect["flag_detect_"+genotype] = np.where(
            (assembFlagMeanDetect["flag_detect_ccs_" + genotype]==1)
            |(assembFlagMeanDetect["flag_detect_shrtRd_rsem_" + genotype]==1)
            |(assembFlagMeanDetect["flag_detect_shrtRd_genomeAln_" + genotype]==1),
            1,
            0
        )
        detectFlagCols.append("flag_detect_"+genotype)

    # Count the number of genotypes detected for each gene and
    #    flag genes detected in at least one genotype but not all (>0, <5)
    assembFlagMeanDetect["num_detect_genotype"] = (
            assembFlagMeanDetect[detectFlagCols].sum(axis=1)
    )
    assembFlagMeanDetect["flag_detect_num_geno_gt0_lt5"] = np.where(
        (assembFlagMeanDetect["num_detect_genotype"] > 0)
        &(assembFlagMeanDetect["num_detect_genotype"] < 5),
        1,
        0
    )

    # Count number of genes with at least one genotype detected, but not all (>0, <5)
#    len(assembFlagMeanDetect[
#        (assembFlagMeanDetect["flag_detect_num_geno_gt0_lt5"]==1)
#    ])
    # 81

    # Count how many of these genes have previously been identified as PAV (flag_pav)
#    len(assembFlagMeanDetect[
#        (assembFlagMeanDetect["flag_detect_num_geno_gt0_lt5"]==1)
#        &(assembFlagMeanDetect["flag_pav"]==1)
#    ])
    # 40

    # Output file of flags and means
    outCols = (
        [
            "gene_id",
            "flag_zea_mays_paralog_hoopes",
            "orthogroup_hoopes",
            "freq_pav_hoopes",
            "flag_pav",
        ]
        + [c for c in assembFlagMeanDetect.columns if "flag_detect" in c]
        + [c for c in assembFlagMeanDetect.columns if "num_detect" in c]
        + [c for c in assembFlagMeanDetect.columns if "analyze" in c]
        + [c for c in assembFlagMeanDetect.columns if "mean" in c]
    )
    assembFlagMeanDetect[outCols].to_csv(args.outFile, index=False)    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
