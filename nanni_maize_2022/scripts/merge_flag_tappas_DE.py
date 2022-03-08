#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge tappAS DEA output files with detection flags")

    # Input data
    parser.add_argument("-d", "--input-directory", dest="inDir", required=True, help="Input directory of maize tappAS output")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file")

    args = parser.parse_args()
    return args

def main():
    # Merge genotype tappAS output files        
    DEAmerge = pd.DataFrame()
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        tempDF = pd.read_csv(
                "{}/{}_tappAS_DEA_Genes.tsv".format(args.inDir,genotype),
                sep="\t"
        )
        tempDF["flag_DE_"+genotype] = np.where(
                tempDF["DEA Result"]=="DE",
                1,
                0
        )
        tempDF = tempDF.rename(columns={
                "#Gene": "gene_id",
                "(1 - Probability)": "DE_pval_"+genotype,
                "Log2FC": "Log2FC_"+genotype,
                "Ambient MeanExpLevel": "mean_TPM_ambient_"+genotype,
                "Ozone MeanExpLevel": "mean_TPM_ozone_"+genotype}
        )
        # Get groups of FC direction
        tempDF[genotype+"_FC_direction"] = np.where(
                tempDF["Log2FC_"+genotype]>0,
                "UP",
                np.where(
                        tempDF["Log2FC_"+genotype]<0,
                        "DOWN",
                        "zero"
                )
        )
        tempDF = tempDF[[
                "gene_id","DE_pval_"+genotype,
                "flag_DE_"+genotype,
                "Log2FC_"+genotype,
                genotype+"_FC_direction",
                "mean_TPM_ambient_"+genotype,
                "mean_TPM_ozone_"+genotype]]
        if genotype == "B73":
            tempMerge = tempDF.copy()
        else:
            tempMerge = pd.merge(
                    DEAmerge,
                    tempDF,
                    how="outer",
                    on="gene_id",
                    indicator="merge_check",
                    validate="1:1"
            )
    #        tempMerge["merge_check"].value_counts()
            del(tempMerge["merge_check"])
        DEAmerge = tempMerge.copy()

    # Make extra flags for GO
    for num in range(1,5):
        DEAmerge["flag_DE_"+str(num)+"_geno"] = np.where(
            DEAmerge[
                [c for c in DEAmerge.columns if "flag_DE_" in c]
            ].sum(axis=1)==num,
            1,
            0
        )
    DEAmerge["flag_DE_all"] = np.where(
            (DEAmerge["flag_DE_B73"]==1)&
            (DEAmerge["flag_DE_C123"]==1)&
            (DEAmerge["flag_DE_Mo17"]==1)&
            (DEAmerge["flag_DE_Hp301"]==1)&
            (DEAmerge["flag_DE_NC338"]==1),
            1,
            0
    )
    DEAmerge["flag_DE_NC338_only"] = np.where(
            (DEAmerge["flag_DE_B73"]==0)&
            (DEAmerge["flag_DE_C123"]==0)&
            (DEAmerge["flag_DE_Mo17"]==0)&
            (DEAmerge["flag_DE_Hp301"]==0)&
            (DEAmerge["flag_DE_NC338"]==1),
            1,
            0
    )
    DEAmerge["flag_DE_B73_only"] = np.where(
            (DEAmerge["flag_DE_B73"]==1)&
            (DEAmerge["flag_DE_C123"]==0)&
            (DEAmerge["flag_DE_Mo17"]==0)&
            (DEAmerge["flag_DE_Hp301"]==0)&
            (DEAmerge["flag_DE_NC338"]==0),
            1,
            0
    )
    DEAmerge["flag_DE_all_noB73"] = np.where(
            (DEAmerge["flag_DE_B73"]==0)&
            (DEAmerge["flag_DE_C123"]==1)&
            (DEAmerge["flag_DE_Mo17"]==1)&
            (DEAmerge["flag_DE_Hp301"]==1)&
            (DEAmerge["flag_DE_NC338"]==1),
            1,
            0
    )

    # Count number of genes DE in at least one genotype
    deDF = DEAmerge[
            (DEAmerge["flag_DE_B73"]==1)|
            (DEAmerge["flag_DE_C123"]==1)|
            (DEAmerge["flag_DE_Mo17"]==1)|
            (DEAmerge["flag_DE_Hp301"]==1)|
            (DEAmerge["flag_DE_NC338"]==1)].copy()
    print("{} gene DE in at least one genotype".format(len(deDF)))

    # Count number of genes with consistent up FC
    cols = [c for c in deDF if "Log2FC" in c]
    upDF = deDF[(deDF[cols] > 0).all(axis=1)]
    downDF = deDF[(deDF[cols] < 0).all(axis=1)]
    print(
        (
            "{0:.2%} ({1} genes total, {2} up, {3} down) have consistent "
            "direction of FC"
        ).format(
                (len(upDF)+len(downDF))/len(deDF),
                len(upDF)+len(downDF),
                len(upDF),len(downDF))
    )

    # Count groups of FC direction
    FCcols = [c for c in deDF.columns if "FC_direction" in c]
    groupDF = deDF.groupby(FCcols)["gene_id"].count().reset_index().sort_values(
            ["gene_id"],ascending=[False]).rename(
                    columns={"gene_id": "num_gene"})

#    groupDF.to_csv(args.outGroup,index=False)

    # Make flags for groups of interest
    DEAmerge["flag_B73Up_restDown"] = np.where(
            ((DEAmerge["flag_DE_B73"]==1)|
                (DEAmerge["flag_DE_C123"]==1)|
                (DEAmerge["flag_DE_Mo17"]==1)|
                (DEAmerge["flag_DE_Hp301"]==1)|
                (DEAmerge["flag_DE_NC338"]==1))&
            (DEAmerge["B73_FC_direction"]=="UP")&
            (DEAmerge["C123_FC_direction"]=="DOWN")&
            (DEAmerge["Hp301_FC_direction"]=="DOWN")&
            (DEAmerge["Mo17_FC_direction"]=="DOWN")&
            (DEAmerge["NC338_FC_direction"]=="DOWN"),
            1,
            0
    )
    DEAmerge["flag_B73Zero_restUp"] = np.where(
            ((DEAmerge["flag_DE_B73"]==1)|
                (DEAmerge["flag_DE_C123"]==1)|
                (DEAmerge["flag_DE_Mo17"]==1)|
                (DEAmerge["flag_DE_Hp301"]==1)|
                (DEAmerge["flag_DE_NC338"]==1))&
            (DEAmerge["B73_FC_direction"]=="zero")&
            (DEAmerge["C123_FC_direction"]=="UP")&
            (DEAmerge["Hp301_FC_direction"]=="UP")&
            (DEAmerge["Mo17_FC_direction"]=="UP")&
            (DEAmerge["NC338_FC_direction"]=="UP"),
            1,
            0
    )
    DEAmerge["flag_B73Down_restUp"] = np.where(
            ((DEAmerge["flag_DE_B73"]==1)|
                (DEAmerge["flag_DE_C123"]==1)|
                (DEAmerge["flag_DE_Mo17"]==1)|
                (DEAmerge["flag_DE_Hp301"]==1)|
                (DEAmerge["flag_DE_NC338"]==1))&
            (DEAmerge["B73_FC_direction"]=="DOWN")&
            (DEAmerge["C123_FC_direction"]=="UP")&
            (DEAmerge["Hp301_FC_direction"]=="UP")&
            (DEAmerge["Mo17_FC_direction"]=="UP")&
            (DEAmerge["NC338_FC_direction"]=="UP"),
            1,
            0
    )
    DEAmerge["flag_NC338Down_restUp"] = np.where(
            ((DEAmerge["flag_DE_B73"]==1)|
                (DEAmerge["flag_DE_C123"]==1)|
                (DEAmerge["flag_DE_Mo17"]==1)|
                (DEAmerge["flag_DE_Hp301"]==1)|
                (DEAmerge["flag_DE_NC338"]==1))&
            (DEAmerge["B73_FC_direction"]=="UP")&
            (DEAmerge["C123_FC_direction"]=="UP")&
            (DEAmerge["Hp301_FC_direction"]=="UP")&
            (DEAmerge["Mo17_FC_direction"]=="UP")&
            (DEAmerge["NC338_FC_direction"]=="DOWN"),
            1,
            0
    )
    print("\n{}\n".format(
        DEAmerge[[c for c in DEAmerge.columns if "flag" in c]].sum().to_string()
    ))
    # Output merged file of flags and tappas results
    DEAmerge.to_csv(args.outFile,index=False)
            
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

