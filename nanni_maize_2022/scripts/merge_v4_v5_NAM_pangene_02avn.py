#!/usr/bin/env python3

# 1) Merge 12604 genes (from supp file 5) with MaizeGDB B73v4 to B73v5 file
# 2) Merge 12604 with NAM pangene list
# 3) Merge with Nature Genetics Mo17-B73 synteny list

import pandas as pd
import numpy as np

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
            "attributes"
        )
    if not gtf["attribute_values"].apply(lambda x: "gene_id" in x).all():
        print(
            "ERROR: gene_id not contained within all exon "
            "attributes"
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

# Get input files
geneDF = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/resubmission_2021/supplement/Supplementary_File_5.csv")

v4tov5 = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/maizeGDB_pan_genes/B73v4_to_B73v5.tsv",
                     sep="\t", names=["B73v4", "B73v5"])

pangene = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/hufford_2021_NAM/pan_gene_matrix_v3_cyverse.csv", low_memory=False)

# Synteny list from Nature Genetics
# Single column of gene ids, using gene as the column name to match AMM in B73v4_Mo17CAU_synteny_NatureGenetics.sas
synListNG = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/Mo17_Table2_gene_list/1.B73_Syntenic_genes/19.B73-sytenic_gene.txt",
                       names=["gene"])
# Add synteny flag
synListNG["flag_Mo17_B73_synteny_NatureGenetics"] = 1


### NAM pangene counts prior to merging ####

# Count total number of pangenes across all 26 NAM lines
len(pangene)
# 103033 pangenes across all 26 NAM lines

# Count total number of pangenes in at least one temperate line
temperate = ["B73", "B97", "HP301", "Il14H", "Ky21", "M162W", "Ms71", "Oh7B", "Oh43", "P39"]
len(pangene[(~pangene[temperate].isna()).any(axis=1)])
# 81097 pangenes in at least one temperate line

# Count tandem duplicates in NAM pangenes (separated by ; under each genotype)
# For all NAM lines
allNAM = [c for c in pangene.columns if c not in ['Pan_gene_ID', 'Pan_gene_id', 'Rep_transcript', 'class', 'Subgenome', 'number_genome_presence']]
pangene["flag_any_tandem_dup"] = pangene[allNAM].apply(
        lambda x: x.str.contains(";"), axis=1).any(axis=1).astype(int)
pangene["flag_any_tandem_dup"].sum()
# 16751 pangenes with at least one tandem duplicate in all NAM lines (~16%)

# For temperate lines
pangene["flag_temperate_tandem_dup"] = pangene[temperate].apply(
        lambda x: x.str.contains(";"), axis=1).any(axis=1).astype(int)
pangene["flag_temperate_tandem_dup"].sum()
# 10729 pangenes with at least one tandem duplicate in any temperate lines (~13% of pangene in temperate)

# Split B73 tandem duplicates
pangeneSplit = split_column_by_sep(pangene, col_name="B73", sep=";")
# Get B73 gene names from transcript names (geneName_T#)
# NOTE: some B73 loci are names with gmap coordiantes (do not change these names)
pangeneSplit["B73_gene_id"] = np.where(
        pangeneSplit["B73"].str.contains("gmap"),
        pangeneSplit["B73"],
        pangeneSplit["B73"].str.split("_").str[0]
)

# Get unique gene_ids with tandem duplicat flags and associated pangene classes
pangeneDF = pangeneSplit.groupby("B73_gene_id").agg({
        "flag_any_tandem_dup": "max",
        "flag_temperate_tandem_dup": "max",
        "class": lambda x: ";".join(x)
    }).reset_index()

### Begin dataset merging ###

# Merge v5 IDs to 12604 by gene_id and v5 IDs
geneV4V5merge = pd.merge(
    geneDF,
    v4tov5,
    how="outer",
    left_on="gene_id",
    right_on="B73v4",
    validate="1:1",
    indicator="merge_check"
)

geneV4V5merge["merge_check"].value_counts()

#right_only    25959
#both          12397
#left_only       207

# 207 (~1.6% of genes) are not in the B73v4 list because they are from ensembl
#       after spot checking some, there are a few that have been dropped from v5

geneV4V5merge["flag_no_assoc_v5_gene"] = np.where(
        geneV4V5merge["merge_check"]=="left_only",
        1,
        0
)
geneV4V5 = geneV4V5merge[geneV4V5merge["merge_check"]!="right_only"].copy().drop(columns=["merge_check"])

# Split the v5 comma separated lists into 1 row each
geneV4V5split = split_column_by_sep(geneV4V5, col_name="B73v5", sep=",")

## Merge with the NAM pangene list
pangeneV4V5merge = pd.merge(
    geneV4V5split,
    pangeneDF,
    how="outer",
    left_on="B73v5",
    right_on="B73_gene_id",
    validate="m:1",
    indicator="merge_check"
)

pangeneV4V5merge["merge_check"].value_counts()

# right_only    39892
# both          13624
# left_only       433
pangeneV4V5merge["flag_no_assoc_pangene"] = np.where(
        pangeneV4V5merge["merge_check"]=="left_only",
        1,
        0
)
pangeneV4V5 = pangeneV4V5merge[pangeneV4V5merge["merge_check"]!="right_only"].copy().drop(columns="merge_check")

# Count the number of unique gene_id Pan_gene_ID pairs
# (if some genes are associated with more than one pangene then there will be more than 12604)
len(pangeneV4V5[["gene_id", "B73_gene_id"]].drop_duplicates())
#14057

# Flag the gene classes
pangeneV4V5["flag_core"] = np.where(pangeneV4V5["class"].fillna("none").apply(lambda x: "Core Gene" in x.split(";")), 1, 0)
pangeneV4V5["flag_near_core"] = np.where(pangeneV4V5["class"].fillna("none").apply(lambda x: "Near-Core Gene" in x.split(";")), 1, 0)
pangeneV4V5["flag_dispensable"] = np.where(pangeneV4V5["class"].fillna("none").apply(lambda x: "Dispensable Gene" in x.split(";")), 1, 0)
pangeneV4V5["flag_private"] = np.where(pangeneV4V5["class"].fillna("none").apply(lambda x: "Private Gene" in x.split(";")), 1, 0)

classFlags = ["flag_core", "flag_near_core", "flag_dispensable", "flag_private"]
dupFlags = ["flag_any_tandem_dup", "flag_temperate_tandem_dup"]

# Get the max of the class flags for each 12604
#   more than one flag means more than one pangene association (multiple v5 gene associations)
pangeneClass = pangeneV4V5.groupby("gene_id")[["flag_no_assoc_v5_gene"] + ["flag_no_assoc_pangene"] + classFlags + dupFlags].max().reset_index()

# Count the class flags - NOTE these are not mutually exclusive due to multiple v5 associated with each v4
pangeneClass[classFlags].sum()
# flag_core           11194
# flag_near_core        782
# flag_dispensable      405
# flag_private           29

# Count the number of the 12604 genes that are not able to be associated with the pangene file
#   they either could not be associated with a v5 ID (mostly ensemble IDs not in NAM annotations)
#   or they are not associated with a pangene in the pandgene file
len(pangeneClass[pangeneClass[classFlags].sum(axis=1)==0])
#433

# Count how many of the 433 are B73-Mo17 syntenic
len(pangeneClass[
    (pangeneClass[classFlags].sum(axis=1)==0)
    & (pangeneClass["gene_id"].isin(synListNG["gene"]))]
)
# 210 of the 433 are B73-Mo17 syntenic

# Count number with an associated pangene
len(pangeneClass[pangeneClass[classFlags].sum(axis=1)>0])
# 12171

# Count the number of the 12604 that are core or near-core
len(pangeneClass[pangeneClass["flag_core"]+pangeneClass["flag_near_core"]>0])
#11936
# (~95% of 12604 or ~97% of 12171 with associated pangene)

# Classify 12604 by the highest level of associated classes
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

# Count the hestest level classes of the 12604
pangeneClass["highest_level_class"].value_counts()

# core           11194
# near-core        742
# none             433
# dispensable      217
# private           18

# Count tandem duplicates in NAM pangenes (separated by ; under each genotype)
# For all NAM lines
pangeneClass["flag_any_tandem_dup"].sum()
# 3680 of the 12171 genes with at least one tandem duplicate in all NAM lines (~30%)
# For temperate lines
pangeneClass["flag_temperate_tandem_dup"].sum()
# 2668 of the 12171 with at least one tandem duplicate in any temperate lines (~22%)

# Get crosstab of pangene class with tandem duplicates
pd.crosstab(pangeneClass["highest_level_class"], pangeneClass["flag_any_tandem_dup"])
# flag_any_tandem_dup   0.0   1.0
# highest_level_class            
# core                 7832  3362
# dispensable           138    79
# near-core             503   239
# private                18     0
pd.crosstab(pangeneClass["highest_level_class"], pangeneClass["flag_temperate_tandem_dup"])
# flag_temperate_tandem_dup   0.0   1.0
# highest_level_class                  
# core                       8751  2443
# dispensable                 162    55
# near-core                   572   170
# private                      18     0

# Merge pangeneClass variables with SFile 5
mergePanSFile5 = pd.merge(
        geneDF,
        pangeneClass,
        how="outer",
        on="gene_id",
        validate="1:1",
        indicator="merge_check"
)
mergePanSFile5["merge_check"].value_counts()
#    both          12604
#    right_only        0
#    left_only         0

# Output new Sfile5 with pangene and tandem duplicate flags
# Rename highest_level_class to be pangene_class
mergePanSFile5.drop(columns=["merge_check"]).rename(
        columns={"highest_level_class": "pangene_class"}).to_csv(
                "~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/genetics_2022/supplement/Supplementary_File_5.csv",
                index=False
)

# Check what are the private genes and which genotypes are they in
privateV4V5 = pangeneV4V5[pangeneV4V5["gene_id"].isin(
        pangeneClass[pangeneClass["highest_level_class"]=="private"]["gene_id"])
]
privateDF = pangeneSplit[
        pangeneSplit["B73_gene_id"].isin(
            privateV4V5["B73_gene_id"]
        )
]
# All are B73 private...DUH used B73 IDs so they will all be B73 private

# Check expression detection/analyzability of the private in our 5 genotypes
privateV4V5[["gene_id", "num_detect_genotype"]].drop_duplicates()["num_detect_genotype"].value_counts().sort_index()
#    4.0     1
#    5.0    17
privateV4V5[["gene_id", "sum_analyze"]].drop_duplicates()["sum_analyze"].value_counts().sort_index()
#    0.0    4
#    1.0    4
#    2.0    6
#    3.0    2
#    4.0    1
#    5.0    1


# Check what are the dispensable genes and which genotypes are they in
dispensableV4V5 = pangeneV4V5[pangeneV4V5["gene_id"].isin(
        pangeneClass[pangeneClass["highest_level_class"]=="dispensable"]["gene_id"])
]
dispensableDF = pangeneSplit[
        pangeneSplit["B73_gene_id"].isin(
            dispensableV4V5["B73_gene_id"]
        )
]

# Count how many dispensable are in HP301
dispensableV4V5[dispensableV4V5["B73_gene_id"].isin(dispensableDF[~dispensableDF["HP301"].isna()]["B73_gene_id"])]["gene_id"].nunique()
# 120
# Count how many dispensable and in HP301 are also in B73-Mo17 syntenic list
dispensableV4V5[
    (dispensableV4V5["B73_gene_id"].isin(dispensableDF[~dispensableDF["HP301"].isna()]["B73_gene_id"]))
    & (dispensableV4V5["gene_id"].isin(synListNG["gene"]))]["gene_id"].nunique()
# 77 of the 120 are also B73-Mo17 syntenic

# Get the genes that are not detected in all 5 genotypes and find NAM classifications
detectLT5 = pangeneV4V5[pangeneV4V5["num_detect_genotype"]<5]
detectLT5["gene_id"].nunique()
# 81
pangeneClass[pangeneClass["gene_id"].isin(detectLT5["gene_id"])][["gene_id", "highest_level_class"]].drop_duplicates()["highest_level_class"].value_counts()
#    none           32
#    core           27
#    dispensable    17
#    near-core       4
#    private         1



# Merge 12604 pangene classes with Nature Genetics Mo17-B73 synteny list
pangeneSynMerge = pd.merge(
        synListNG,
        pangeneClass,
        how="outer",
        left_on = "gene",
        right_on = "gene_id",
        validate = "1:1",
        indicator = "merge_check"
)
pangeneSynMerge["merge_check"].value_counts()
#    left_only     21803
#    both          11878
#    right_only      726

# Drop the left_only, or genes only in the synteny list
pangeneSyn = pangeneSynMerge[pangeneSynMerge["merge_check"]!="left_only"].drop(columns=["merge_check"])

# Get counts of pangene classes in genes that are not in syntenic list
len(pangeneSyn[pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1])
# 726
pangeneSyn[pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1]["highest_level_class"].value_counts()
#    core           361
#    none           223
#    dispensable    100
#    near-core       29
#    private         13

# Which of the near-core are in Hp301
nearCoreNoSyn = pangeneSyn[
        (pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1)
        & (pangeneSyn["highest_level_class"]=="near-core")
].copy()
dispensableDF = pangeneSplit[
        pangeneSplit["B73_gene_id"].isin(
            dispensableV4V5["B73_gene_id"]
        )
]

nearCoreNoSynFULL = pangeneV4V5[pangeneV4V5["gene_id"].isin(nearCoreNoSyn["gene_id"])]
nearCoreNoSynFULL2 = pangeneSplit[pangeneSplit["B73_gene_id"].isin(nearCoreNoSynFULL["B73_gene_id"])]
nearCoreNoSynFULL[nearCoreNoSynFULL["B73_gene_id"].isin(nearCoreNoSynFULL2[~nearCoreNoSynFULL2["HP301"].isna()]["B73_gene_id"])]["gene_id"].nunique()
# 29
# All near-core not in the Mo17-B73 are in Hp301

# Which of the dispensable are in Hp301
dispNoSyn = pangeneSyn[
        (pangeneSyn["flag_Mo17_B73_synteny_NatureGenetics"]!=1)
        & (pangeneSyn["highest_level_class"]=="dispensable")
].copy()
dispNoSynFULL = pangeneV4V5[pangeneV4V5["gene_id"].isin(dispNoSyn["gene_id"])]
dispNoSynFULL2 = pangeneSplit[pangeneSplit["B73_gene_id"].isin(dispNoSynFULL["B73_gene_id"])]
dispNoSynFULL[dispNoSynFULL["B73_gene_id"].isin(dispNoSynFULL2[~dispNoSynFULL2["HP301"].isna()]["B73_gene_id"])]["gene_id"].nunique()
# 43
# 43 of the 100 dispensable not in the Mo17-B73 are in Hp301



# What is the expression of the near-core in other genotypes (detection flags)
nearCoreNoSynFULL[["gene_id", "num_detect_genotype"]].drop_duplicates()["num_detect_genotype"].value_counts()
# 5.0    29
# All 29 near-core not in syntenic genes are detected in all 5 genotypes in long/short read data
nearCoreNoSynFULL[["gene_id", "sum_analyze"]].drop_duplicates()["sum_analyze"].value_counts().sort_index()
#    0.0     1
#    1.0     2
#    2.0     1
#    5.0    25
analyzeFlags = [c for c in nearCoreNoSynFULL.columns if "flag_analyze" in c]
nearCoreNoSynFULL3 = nearCoreNoSynFULL[["gene_id", "sum_analyze"] + analyzeFlags].drop_duplicates()
print(nearCoreNoSynFULL3[nearCoreNoSynFULL3["sum_analyze"]<5].to_string(index=False))
#        gene_id  sum_analyze  flag_analyze_B73  flag_analyze_C123  flag_analyze_Hp301  flag_analyze_Mo17  flag_analyze_NC338
# Zm00001d000223          1.0               0.0                1.0                 0.0                0.0                 0.0
# Zm00001d023373          2.0               0.0                1.0                 0.0                0.0                 1.0
# Zm00001d023796          0.0               0.0                0.0                 0.0                0.0                 0.0
# Zm00001d042625          1.0               0.0                0.0                 0.0                1.0                 0.0

# 1 gene not analyzable in any genotype
# 2 genes analyzable in 1 genotype (C123 for one and Mo17 the other)
# 1 gene analyzable in 2 genotypes (C123 and NC338)


# What is the expression of the dispensable in other genotypes (detection flags)
dispNoSynFULL[["gene_id", "num_detect_genotype"]].drop_duplicates()["num_detect_genotype"].value_counts().sort_index()
#    1.0     1
#    2.0     1
#    3.0     4
#    4.0     6
#    5.0    88
# Nearly all (88 of 100 dispensable not in syntenic genes are detected in all 5 genotypes in long/short read data
dispNoSynFULL[["gene_id", "sum_analyze"]].drop_duplicates()["sum_analyze"].value_counts().sort_index()
#    0.0    14
#    1.0    27
#    2.0    31
#    3.0    14
#    4.0     6
#    5.0     8


# Are they paralogs of something else



## How many of the 789 (788 in same direction) DE in all non-B73 genotypes are tandem duplicates
deDF = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/resubmission_2021/supplement/Supplementary_File_4.csv", low_memory=False)
deCols = [c for c in deDF.columns if "flag_analyze_DE" in c]
nonB73de = deDF[(deDF["flag_analyze_DE_B73"]==0)&(deDF[deCols].sum(axis=1)==4)]
deFULL = pangeneV4V5[pangeneV4V5["gene_id"].isin(nonB73de["gene_id"])]
deTandemDups = deFULL.groupby("gene_id")["flag_any_tandem_dup"].max().reset_index()
deTandemDups["flag_any_tandem_dup"].sum()
# 267

# Count how many of the 788 that are the same direction
logCols = [c for c in deDF.columns if "Log2FC" in c and "B73" not in c]
nonB73deSame = nonB73de[(nonB73de[logCols]>0).all(axis=1)|(nonB73de[logCols]<0).all(axis=1)]
deSameFULL = pangeneV4V5[pangeneV4V5["gene_id"].isin(nonB73deSame["gene_id"])]
deSameTandemDups = deSameFULL.groupby("gene_id")["flag_any_tandem_dup"].max().reset_index()
deSameTandemDups["flag_any_tandem_dup"].sum()
# 267

# How many of the 12604 are orthologs of sorghum
# !!! NOTE: Need to get this comparison for the orthologs in the newest 2021 paper

# Merge 12604 with hoopes paralog/ortholog file
orthoDF = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/Hoopes_2018/Zea_mays_Hoopes_2018_rice_sorg_arabidop_amborel_ortho_para_flag.csv")
geneOrthoMerge = pd.merge(
        orthoDF,
        geneDF,
        how = "outer",
        on = "gene_id",
        validate = "1:1",
        indicator = "merge_check"
)
geneOrthoMerge["merge_check"].value_counts()
#    left_only     16612
#    both           7846
#    right_only     4758

# Get only those in 12604 (not left_only)
geneOrtho = geneOrthoMerge[geneOrthoMerge["merge_check"]!="left_only"].copy()

# Count number of 12604 with each ortho flag type
geneOrtho[[c for c in orthoDF.columns if "flag" in c]].sum()
#    flag_Zea_mays_paralog              4473.0
#    flag_ortho_Amborella_trichopoda    6057.0
#    flag_ortho_Arabidopsis_thaliana    6182.0
#    flag_ortho_Oryza_sativa            6914.0
#    flag_ortho_Sorghum_bicolor         7224.0


# Merge with all B73 v4 genes
allGene = get_gene_from_gtf("~/mclab/SHARE/McIntyre_Lab/useful_maize_info/Zea_mays.B73_RefGen_v4.41.gtf")
allGeneOrthoMerge = pd.merge(
        geneOrthoMerge,
        allGene,
        how = "outer",
        on = "gene_id",
        validate = "1:1",
        indicator = "merge_check2"
)
# Check that all 12604 are present (no left_only)
allGeneOrthoMerge["merge_check2"].value_counts()
#    both          29216
#    right_only    17214
#    left_only         0

# Flag 12604
allGeneOrthoMerge["flag_in_12604"] = np.where(
        allGeneOrthoMerge["gene_id"].isin(geneDF["gene_id"]),
        1,
        0
)

# Get crosstab of 12604 and sorghum orthologs in all B73 v4 genes
pd.crosstab(allGeneOrthoMerge["flag_in_12604"], allGeneOrthoMerge["flag_ortho_Sorghum_bicolor"].fillna(0))
#    flag_ortho_Sorghum_bicolor    0.0    1.0
#    flag_in_12604                           
#    0                           23058  10768
#    1                            5380   7224

# ~61% of 12604 have a sorghum ortholog
# ~32% of non-12604 reference genes have a sorghum ortholog
# ~39% of all reference genes (17992 / 46430) have a sorghum ortholog
