#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sqlite3
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in merged flags of detection above input in all sample types within a species to make counts")

    # Input data
    parser.add_argument("-i", "--input-flags", dest="inFlags", required=True, help="Input CSV of merged flag file, has featureID and a column of each sampleType (antiboy_species_genotype_sex_treatment) ")
    parser.add_argument("-f", "--featureTypes", dest="inFeature", required=True, help="Concatenated CSV file connecting featureID to featureType")
    parser.add_argument("-s", "--species", dest="inSpecies", required=True, help="Species name")
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temp SQL database file to reduce memory required", default=":memory:")
    parser.add_argument("-p", "--feature-path",dest="featPath",required=False, help="Directory where feature files with associated FBgn are located (*_xcrpt_annotation.csv)")

    # Output data
    parser.add_argument("-o", "--output", dest="outPrefix", required=True, help="Prefix for output files")
    parser.add_argument("-l", "--listDir", dest="outDir", required=False, help="Output directory for lists of genes in each flag category, directory must be made prior to running")

    args = parser.parse_args()
    
    # Dependent arguments
    if args.featPath and not args.outDir:
        parser.error("--feature-path argument requires --outdir argument")
    elif args.outDir and not args.featPath:
        parser.error("--outdir argument requires --feature-path argument")
    return args
    
def main():
    
    # Get stacked flag file and featureID to featureType file
    stackedDF = pd.read_csv(args.inFlags)
    featureDF = pd.read_csv(args.inFeature)

    # Connect to SQL database
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()

    # Merge feature type into flag file
    stackedDF.to_sql("flags", con, if_exists="replace")
    featureDF.to_sql("features", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in2.featureType, in1.* "
                "FROM flags in1 LEFT JOIN features in2 "
                "ON in1.featureID = in2.featureID ;")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])

    # Flag female and male on
    f4Col = [col for col in mergeDF.columns if ("_f_" in col)&("K4_" in col)]
    f27Col = [col for col in mergeDF.columns if ("_f_" in col)&("K27_" in col)]
    m4Col = [col for col in mergeDF.columns if ("_m_" in col)&("K4_" in col)]
    m27Col = [col for col in mergeDF.columns if ("_m_" in col)&("K27_" in col)]
    mergeDF['flag_f_K4_on'] = np.where(mergeDF[f4Col].sum(axis=1)==2,1,0)
    mergeDF['flag_f_K27_on'] = np.where(mergeDF[f27Col].sum(axis=1)==2,1,0)
    mergeDF['flag_m_K4_on'] = np.where(mergeDF[m4Col].sum(axis=1)==2,1,0)
    mergeDF['flag_m_K27_on'] = np.where(mergeDF[m27Col].sum(axis=1)==2,1,0)
    mergeDF.to_csv(args.outPrefix+"_flag.csv",index=False)

    # Get TSS300bpWindow X vs. autosome crosstabs for each flag
    # flag_x = 1 when on X, 0.5 when on chr 4, -1 when scaffold, and 0 when autosome
    TSSDF = mergeDF[mergeDF['featureType']=="TSS300bpWindow"]
    TSSDF['flag_x'] = np.where((TSSDF['featureID'].str.contains("X"))|(TSSDF['featureID'].str.contains("Scf_X")),"1",
                    np.where((TSSDF['featureID'].str.contains("TSS300bpWindow_4_"))|(TSSDF['featureID'].str.contains("TSS300bpWindow_Scf_4_")),"0.5",
                    np.where((TSSDF['featureID'].str.contains("Scaffold"))|(TSSDF['featureID'].str.contains("NODE")),"-1","0")))
    TSSDF.to_csv(args.outPrefix+"_TSS300bpWindow_x_auto.csv", index=False)
    crossOut = open(args.outPrefix+"_TSS300bpWindow_x_auto_crosstabs.txt",'w')
    crossOut.write("## Flag information:\n\tflag_x==1 when on X\n\tflag_x==0.5 when on chr 4\n\tflag_x==0 when on autosome\n\tflag_x==-1 when on scaffold/contig\n\n")
    crossOut.write(pd.crosstab(TSSDF['flag_f_K4_on'],TSSDF['flag_x']).to_string()+"\n\n")
    crossOut.write(pd.crosstab(TSSDF['flag_m_K4_on'],TSSDF['flag_x']).to_string()+"\n\n")
    crossOut.write(pd.crosstab(TSSDF['flag_f_K27_on'],TSSDF['flag_x']).to_string()+"\n\n")
    crossOut.write(pd.crosstab(TSSDF['flag_m_K27_on'],TSSDF['flag_x']).to_string()+"\n\n")
    crossOut.close()

    # Output flag counts
    countDF = mergeDF.groupby('featureType')[['flag_f_K27_on','flag_f_K4_on','flag_m_K27_on','flag_m_K4_on']].sum().reset_index()
    countDF.to_csv(args.outPrefix+"_flag_counts.csv",index=False)

    # Output merged flag string counts (crosstab and individual counts)
    t = open(args.outPrefix+"_fK4_fK27_mK4_mK27_crosstab.csv",'w')
    flags = mergeDF.set_index('featureID').reindex(mergeDF.filter(regex="(flag_|featureType)").columns,axis=1)
    flags['sum_type_flags'] = flags.astype(str).values.sum(axis=1)
    flags['sum_flags'] = flags[[col for col in flags.columns if "flag_" in col]].astype(str).values.sum(axis=1)
    t.write(pd.crosstab(flags['featureType'],flags['sum_flags']).to_string()+"\n")
    t.close()
    catFlags = flags['sum_type_flags'].value_counts().sort_index().rename('featureType_fK4_fK27_mK4_mK27')
    catFlags.reset_index().to_csv(args.outPrefix+"_fK4_fK27_mK4_mK27_counts.csv",index=False)

    # Check if outDir was given
#    if args.outDir != "":
#        # Get different groups of male or female , K4 or K27 and print featureIDs
#    K27_f_only_noK4
#    K27_m_only_noK4
#    K4_f_only_noK27
#    K4_f_no_K4_m
#    K4_m_only_noK27
#    K4_m_no_K4_f
#    K4_both_sexes_noK27
#    K4_both_sexes
#    K27_both_sexes_noK4
#    K4_both_sexes
#    f_both_antibodies_noM
#    m_both_antibodies_noF
#        # Get different groups of male or female , K4 or K27 and print featureIDs per featureType
#        uniqTypes = flags['featureType'].unique()
#        uniqTypes = uniqTypes[uniqTypes != "intergenic"]
#        for featureType in uniqTypes:
#            # Get featureType annotation file
#            featTOgene = pd.read_csv("{}/{}_{}_xcrpt_annotation.csv".format(args.featPath,args.inSpecies,featureType))

#            # Fix fragment fusion and intron column headers and split the FBgn values
#            # Also output files of all FBgn in fragments and introns
#            if "fragment" in group:
#                tempDF = featTOgene.rename(columns={'gene_id':'FBgn','fragment_id':'featureID'})
#                featTOgene = pd.DataFrame(tempDF['FBgn'].str.split('|').tolist(),index=tempDF['featureID']).stack().reset_index([0,'featureID'])
#                featTOgene.columns = ['featureID','FBgn']
#                featTOgene['FBgn'].drop_duplicates().to_csv("{}/{}_fragment_all.txt".format(args.outDir,args.inSpecies), index=False)
#            elif "fusion" in group:
#                tempDF = featTOgene.rename(columns={'gene_id':'FBgn','fusion_id':'featureID'})
#                featTOgene = pd.DataFrame(tempDF['FBgn'].str.split('|').tolist(),index=tempDF['featureID']).stack().reset_index([0,'featureID'])
#                featTOgene.columns = ['featureID','FBgn']
#                featTOgene['FBgn'].drop_duplicates().to_csv("{}/{}_fusion_all.txt".format(args.outDir,args.inSpecies), index=False)
#            elif "intron" in group:
#                tempDF = featTOgene.rename(columns={'gene_id':'FBgn','intron_id':'featureID'})
#                featTOgene = pd.DataFrame(tempDF['FBgn'].str.split('|').tolist(),index=tempDF['featureID']).stack().reset_index([0,'featureID'])
#                featTOgene.columns = ['featureID','FBgn']
#               featTOgene['FBgn'].drop_duplicates().to_csv("{}/{}_intron_all.txt".format(args.outDir,args.inSpecies), index=False)

#            # Merge featureIDs with table of associated FBgn and print out FBgn list
#                featTOgene.to_sql("features", con, if_exists="replace")
#                pd.DataFrame(flags[flags['sum_type_flags']==group].index).to_sql("groups", con, if_exists="replace")
#                cur.execute("DROP TABLE IF EXISTS genes; ")
#                cur.execute("CREATE TABLE genes AS SELECT in1.featureID, in2.FBgn "
#                            "FROM groups in1 LEFT JOIN features in2 "
#                            "ON in1.featureID = in2.featureID ;")
#                geneDF = pd.read_sql("SELECT * FROM genes", con)
#                filename = "{}/{}_{}_{}.txt".format(
#                        args.outDir,args.inSpecies,group.strip("10"),groupType)
#                geneDF['FBgn'].drop_duplicates().to_csv(filename,index=False)
#                othername = "{}/{}_{}_{}_featureID.csv".format(
#                        args.outDir,args.inSpecies,group.strip("10"),groupType)
#                geneDF[['featureID','FBgn']].drop_duplicates().to_csv(othername,index=False)
#                # Print counts
#                print("\tGroup: {} {}\n\t# featureID: {}\n\t# FBgn: {}".format(
#                      group.strip("10"),groupType,
#                      len(flags[flags['sum_type_flags']==group].index),
#                      len(geneDF['FBgn'].drop_duplicates())))


#        for group in flags['sum_type_flags'].unique():
#            if "intergenic" in group:
#                continue
#            else:
#                if "1000" in group:
#                    groupType = "K4_f_only"
#                elif "0100" in group:
#                    groupType = "K27_f_only"
#                elif "0010" in group:
#                    groupType = "K4_m_only"
#                elif "0001" in group:
#                    groupType = "K27_m_only"
#                elif "0011" in group:
#                    groupType = "m_both_antibody"
#                elif "1100" in group:
#                    groupType = "f_both_antibody"
#                elif "1010" in group:
#                    groupType = "K4_both_sexes"
#                elif "0101" in group:
#                    groupType = "K27_both_sexes"
#                else:
#                    continue             
   
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
