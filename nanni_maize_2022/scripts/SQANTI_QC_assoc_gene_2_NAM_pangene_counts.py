#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # Input data
    parser.add_argument(
        "-c",
        "--classification",
        dest="inClass",
        required=True,
        help="SQANTI QC classification file (combined across samples with sample and genotype columns)."
    )
    parser.add_argument(
        "-v",
        "--v4-to-v5",
        dest="inV",
        required=True,
        help="MaizeGDB B73 v4 to v5 conversion file."
    )
    parser.add_argument(
        "-p",
        "--NAM-pangene",
        dest="inP",
        required=True,
        help="NAM pangene file from Hufford 2021."
    )
#    parser.add_argument(
#        "-s",
#        "--synteny",
#        dest="inS",
#        required=True,
#        help="Synteny list from Nature Genetics."
#    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output file for counts."
    )

    args = parser.parse_args()
    return args

def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
    # Split variable by some character like '|' or ',' and keep all other values the same
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def main():
    # 1) Count unique number of annotated gene and novel loci mapped by each sample
    # 2) Make wide file of cluster count for each annotated gene across samples
    # 3) Merge cluster associated annotated genes with MaizeGDB B73v4 to B73v5 file
    # 4) Merge cluster associated annotated genes with NAM pangene list and count
    # 5) Merge cluster associated annotated genes with Nature Genetics Mo17-B73 synteny list and count

    # Open output file
    outFile = open(args.outFile, "w")

    # Get input files
    classDF = pd.read_csv(args.inClass, sep="\t", low_memory=False)
    v4tov5 = pd.read_csv(args.inV, sep="\t", names=["B73v4", "B73v5"])
#    v4tov5 = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/maizeGDB_pan_genes/B73v4_to_B73v5.tsv",
#                         sep="\t", names=["B73v4", "B73v5"])
    pangene = pd.read_csv(args.inP, low_memory=False)
#    pangene = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/hufford_2021_NAM/pan_gene_matrix_v3_cyverse.csv", low_memory=False)
    
    # Synteny list from Nature Genetics
    # Single column of gene ids, using gene as the column name to match AMM in B73v4_Mo17CAU_synteny_NatureGenetics.sas
#    synListNG = pd.read_csv(args.inS, names=["gene"])
#    synListNG = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/Mo17_Table2_gene_list/1.B73_Syntenic_genes/19.B73-sytenic_gene.txt",
#                           names=["gene"])
    # Add synteny flag
#    synListNG["flag_Mo17_B73_synteny_NatureGenetics"] = 1

    # Add flag for pangene presence in all 10 temperate genotypes (including B73)
    temperate = ["B73", "B97", "HP301", "Il14H", "Ky21", "M162W", "Ms71", "Oh7B", "Oh43", "P39"]
    pangene["flag_in_all_temperate"] = np.where(
            ~pangene[temperate].isna().any(axis=1),
            1,
            0
    )
    outFile.write(
        "Count of NAM pangene classifications to those contained in all temperate "
        "({})\n{}\n\n".format(
                temperate,
                pd.crosstab(pangene["flag_in_all_temperate"],pangene["class"]).to_string()))

    # Get number of clusters per loci for each sample
    sampleDF = classDF.groupby(
            ["sample","genotype","associated_gene"]
        )["isoform"].count().reset_index().rename(
                columns={"isoform": "num_cluster"})

    # 1) Count unique number of annotated gene and novel loci mapped by each sample
    # NOTE: annotated genes in b73 do not contain "_" so only novel loci will have this
    annotDF = sampleDF[~sampleDF["associated_gene"].str.contains("_")]
    annotCount = annotDF.groupby("sample")[
        "associated_gene"].nunique().reset_index().rename(
            columns={"associated_gene": "num_annotated_gene"})
    annotCountSamp = annotDF.groupby("associated_gene")["sample"].count().reset_index()
    annotCountGenotype = annotDF.groupby("associated_gene")["genotype"].nunique().reset_index()
    outFile.write("{}\n\n".format(annotCount.to_string(index=False)))
    outFile.write("{} unique annotated genes across all samples (union)\n\n".format(
            annotDF["associated_gene"].nunique()))
    outFile.write("{} unique annotated genes across all samples (intersection)\n\n".format(
            len(annotCountSamp[annotCountSamp["sample"]==len(annotCount)])))
    outFile.write("{} unique annotated genes across all genotypes (intersection)\n\n".format(
            len(annotCountGenotype[annotCountGenotype["genotype"]==annotDF["genotype"].nunique()])))
    novelCount = sampleDF[
            sampleDF["associated_gene"].str.contains("_")
        ].groupby("sample")["associated_gene"].nunique().reset_index().rename(
            columns={"associated_gene": "num_novel_loci"})
    outFile.write("{}\n\n".format(novelCount.to_string(index=False)))
    novelNoASFusionCount = sampleDF[
            (sampleDF["associated_gene"].str.contains("_"))
            & (sampleDF["associated_gene"].str.contains("novel"))
            & (~sampleDF["associated_gene"].str.endswith("_AS"))
        ].groupby("sample")["associated_gene"].nunique().reset_index().rename(
            columns={"associated_gene": "num_novel_loci_no_AS_or_fusion"})
    outFile.write("{}\n\n".format(novelNoASFusionCount.to_string(index=False)))


    # 2) Make wide file of cluster count for each annotated gene across samples
    wideAnnotDF = sampleDF[
            ~sampleDF["associated_gene"].str.contains("_")
        ].pivot_table(
            index=["associated_gene"],
            columns="sample",
            values="num_cluster").reset_index().fillna(0)


    # 3) Merge cluster associated annotated genes with MaizeGDB B73v4 to B73v5 file
    geneV4V5merge = pd.merge(
        wideAnnotDF,
        v4tov5,
        how="outer",
        left_on="associated_gene",
        right_on="B73v4",
        validate="1:1",
        indicator="merge_check"
    )
    outFile.write(
        "After merge with V4-to-V5 (left-only indicate # missing from V5)\n{}\n\n".format(
                geneV4V5merge["merge_check"].value_counts().to_string()))
    
    geneV4V5 = geneV4V5merge[geneV4V5merge["merge_check"]!="right_only"].copy().drop(columns=["merge_check"])
    
    # Split the v5 comma separated lists into 1 row each
    geneV4V5split = split_column_by_sep(geneV4V5, col_name="B73v5", sep=",")


    # 4) Merge cluster associated annotated genes with NAM pangene list and count

    ## Merge with the NAM pangene list
    pangeneV4V5merge = pd.merge(
        geneV4V5split,
        pangene,
        how="outer",
        left_on="B73v5",
        right_on="Pan_gene_id",
        validate="m:1",
        indicator="merge_check"
    )
    outFile.write(
        "After merge with NAM pangenes (left-only indicate # missing from NAM pangenes)\n{}\n\n".format(
                pangeneV4V5merge["merge_check"].value_counts().to_string()))
    
    pangeneV4V5 = pangeneV4V5merge[pangeneV4V5merge["merge_check"]!="right_only"].copy().drop(columns="merge_check")
    
    # Count the number of unique gene_id Pan_gene_ID pairs
    # (some genes are associated with more than one pangene so there can be duplicates)
    outFile.write(
        "{} unique gene-pangene pairs (some genes are associated with more "
        "than one pangene so there can be duplicates)\n\n".format(
            len(pangeneV4V5[["associated_gene", "Pan_gene_ID"]].drop_duplicates())))
    
    # Flag the gene classes
    pangeneV4V5["flag_core"] = np.where(pangeneV4V5["class"]=="Core Gene", 1, 0)
    pangeneV4V5["flag_near_core"] = np.where(pangeneV4V5["class"]=="Near-Core Gene", 1, 0)
    pangeneV4V5["flag_dispensable"] = np.where(pangeneV4V5["class"]=="Dispensable Gene", 1, 0)
    pangeneV4V5["flag_private"] = np.where(pangeneV4V5["class"]=="Private Gene", 1, 0)
    
    classFlags = ["flag_core", "flag_near_core", "flag_dispensable", "flag_private"]
    
    # Get the max of the class flags for each associated genes
    #   more than one flag means more than one pangene association (multiple v5 gene associations)
    pangeneClass = pangeneV4V5.groupby("associated_gene")[
            classFlags + ["flag_in_all_temperate"]].max().reset_index()
    
    # Count the class flags - NOTE these are not mutually exclusive due to multiple v5 associated with each v4
    outFile.write(
        "Counts of pangene classifications (NOTE: Not mutually exclusive due "
        "to multiple v5 associations with each v4 gene)\n{}\n\n".format(
                pangeneClass[classFlags + ["flag_in_all_temperate"]].sum().to_string()))
    
    # Count the number of the associated genes that are not able to be associated with the pangene file
    #   they either could not be associated with a v5 ID (mostly ensemble IDs not in NAM annotations)
    #   or they are not associated with a pangene in the pandgene file
    outFile.write(
        "{} genes could not be associated with a pangene due to either no v5 "
        "ID or no associated pangene in the NAM pangene list.".format(
                len(pangeneClass[pangeneClass[classFlags].sum(axis=1)==0])))

    # Classify associated genes by the highest level of associated classes
    #   core > near-core > dispensible > private
    classConditions = [
            pangeneClass["flag_core"]==1,
            pangeneClass["flag_near_core"]==1,
            pangeneClass["flag_dispensable"]==1,
            pangeneClass["flag_private"]==1,
    ]
    classChoices = [
            "core",
            "near-core",
            "dispensable",
            "private"
    ]
    pangeneClass["highest_level_class"] = np.select(classConditions, classChoices, "none")

    # Count the highest level classes of the associated genes
    outFile.write(
        "Counts of highest pangene classifications for each associated gene:\n{}\n\n".format(
                pangeneClass["highest_level_class"].value_counts().to_string()))

    # Count the number/proportion of core, near-core, or core/near-core
    outFile.write(
        "{0} associated genes ({1:.2%}) are core\n{2} associated genes "
        "({3:.2%}) are near-core\n{4} associated genes ({5:.2%}) are core or "
        "near-core\n\n".format(
            len(pangeneClass[pangeneClass["highest_level_class"]=="core"]),
            len(pangeneClass[pangeneClass["highest_level_class"]=="core"])/len(pangeneClass),
            len(pangeneClass[pangeneClass["highest_level_class"]=="near-core"]),
            len(pangeneClass[pangeneClass["highest_level_class"]=="near-core"])/len(pangeneClass),
            len(pangeneClass[pangeneClass["highest_level_class"].isin(["core","near-core"])]),
            len(pangeneClass[pangeneClass["highest_level_class"].isin(["core","near-core"])])/len(pangeneClass)))

    # Count the number of genes in all 10 temperate NAM genotypes
    outFile.write(
            "{0} associated genes ({1:.2%}) are in all 10 temperate NAM genotypes".format(
                    len(pangeneClass[pangeneClass["flag_in_all_temperate"]==1]),
                    len(pangeneClass[pangeneClass["flag_in_all_temperate"]==1])/len(pangeneClass)))

    # Check the B73 private genes for how many of our samples found evidence
    privateV4V5 = pangeneV4V5[pangeneV4V5["associated_gene"].isin(
            pangeneClass[pangeneClass["highest_level_class"]=="private"]["associated_gene"])
    ]
    privateSamples = privateV4V5.groupby("associated_gene")[annotCount["sample"]].max()
    privateSamplesGroup = (privateSamples > 0).astype(int).reset_index().groupby(
            list(annotCount["sample"])
        )["associated_gene"].count().reset_index()
    outFile.write(
        "For the {} genes private to B73, they were in the following samples:\n{}\n\n".format(
                len(privateSamples),
                privateSamplesGroup.to_string(index=False)))
    outFile.write("Genes private to B73 and only found in B73 samples:\n{}\n\n".format(
        privateSamples[
                (privateSamples[[c for c in annotCount["sample"] if "b73" in c]].sum(axis=1)>0)
                & (privateSamples[[c for c in annotCount["sample"] if "b73" not in c]].sum(axis=1)==0)
            ].reset_index()[
            [c for c in annotCount["sample"] if "b73" in c] + ["associated_gene"]].to_string(index=False)))



    # 5) Merge cluster associated annotated genes with Nature Genetics Mo17-B73 synteny list and count
        
    # Merge annotated gene pangene classes with Nature Genetics Mo17-B73 synteny list
#    pangeneSynMerge = pd.merge(
#            synListNG,
#            pangeneClass,
#            how="outer",
#            left_on = "gene",
#            right_on = "associated_gene",
#            validate = "1:1",
#            indicator = "merge_check"
#    )
#    
#    # Drop the left_only, or genes only in the synteny list
#    pangeneSyn = pangeneSynMerge[pangeneSynMerge["merge_check"]!="left_only"].drop(columns=["merge_check"])
#    
#    # Get counts of pangene classes in genes that are not in syntenic list
#    outFile.write("{} ".format(
#            len(pangeneSyn[pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1])))
#    # 726
#    pangeneSyn[pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1]["highest_level_class"].value_counts()
#    #    core           335
#    #    none           255
#    #    dispensable     95
#    #    near-core       28
#    #    private         13
#    
#    # Which of the near-core are in Hp301
#    nearCoreNoSyn = pangeneSyn[
#            (pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1)
#            & (pangeneSyn["highest_level_class"]=="near-core")
#    ].copy()
#    nearCoreNoSynFULL = pangeneV4V5[pangeneV4V5["gene_id"].isin(nearCoreNoSyn["gene_id"])]
#    nearCoreNoSynFULL[~nearCoreNoSynFULL["HP301"].isna()]["gene_id"].nunique()
#    
#    # Which of the dispensable are in Hp301
#    dispNoSyn = pangeneSyn[
#            (pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1)
#            & (pangeneSyn["highest_level_class"]=="dispensable")
#    ].copy()
#    dispNoSynFULL = pangeneV4V5[pangeneV4V5["gene_id"].isin(dispNoSyn["gene_id"])]
#    dispNoSynFULL[~dispNoSynFULL["HP301"].isna()]["gene_id"].nunique()


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

