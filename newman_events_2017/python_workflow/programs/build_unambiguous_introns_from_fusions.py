#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-09
# NAME: build_unambiguous_introns_from_fusions.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: Using fusions (exonic regions) create an annotation file of all unambiguous intronic regions
# within a gene
#
# REQUIRED PACKAGES: argparse  (tested with v1.1)
#                    pandas    (tested with v0.19.2)
#                    logging   (tested with v0.5.1.2)
#
#######################################################################################################################

# Packages

import argparse
import logging
import pandas as pd

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import fusion annotation files and create annotations for"
                                                 "unambiguous intronic sequences")
    parser.add_argument("--input-fusion-file", dest="inFusInfo", action='store', required=True,
                        help="Input fusion info TSV file")
    parser.add_argument("--input-fusion-bed", dest="inFusBED", action='store', required=True,
                        help="Input fusion BED file")
    parser.add_argument("--outCSV", dest="outCSV", action='store', required=True, help="Output CSV file name")
    parser.add_argument("--outBED", dest="outBED", action='store', required=True, help="Output BED file name")
    args = parser.parse_args()
    return args

def main():
    # Import fusion annotations
    logger.info("Importing input annotations")
    fusInfoDF = pd.read_csv(args.inFusInfo, sep="\t", skiprows=1,
                            names=['fusion_id','exon_id','gene_id','flag_multigene'],
                            usecols=['fusion_id','gene_id','flag_multigene'])
    fusBedDF = pd.read_csv(args.inFusBED, sep="\t", header=None,
                           names=['chr','fusion_start','fusion_stop','fusion_id','score','strand'],
                            usecols=['chr','fusion_start','fusion_stop','fusion_id'])
    # Sort fusion DFs and merge
    fusInfoDF = fusInfoDF.sort_values(by=['fusion_id'])
    fusInfoDF = fusInfoDF.drop_duplicates(keep='first')
    fusBedDF = fusBedDF.sort_values(by=['fusion_id'])
    fusInfoDF2 = pd.merge(fusBedDF,fusInfoDF, how='inner', on=['fusion_id'])

    # Drop genes with fewer than 2 fusions
    fusPerGene=fusInfoDF2['gene_id'].value_counts().to_frame()
    multiFusPerGene=fusPerGene[fusPerGene.gene_id > 1]
    genesKeep = multiFusPerGene.index.tolist()

    # Set up my lists for intron start, intron stop, ID, chromosome, geneID, and multigeneFlags
    intronIDList=[]
    intronStartList=[]
    intronStopList=[]
    intronSizeList=[]
    intronChrList=[]
    intronGeneList=[]
    intronFus1List=[]
    intronFus2List = []
    intronFus1Multigene=[]
    intronFus2Multigene=[]
    # Iterate over genes
    for gn in range(0,len(genesKeep)):
        geneID=genesKeep[gn]
        # Set up lists
        chrList=[]
        fusStartList=[]
        fusStopList=[]
        fusIdList=[]
        flagMultiGeneList=[]
        # Subset fusion dataframe for fusions from gene
        subsetFusDF = fusInfoDF2[fusInfoDF2['gene_id'] == geneID]
        subsetFusDF = subsetFusDF.sort_values(by=['fusion_start','fusion_stop'])
        chrList=subsetFusDF['chr'].tolist()
        fusStartList=subsetFusDF['fusion_start'].tolist()
        fusStopList=subsetFusDF['fusion_stop'].tolist()
        fusIdList=subsetFusDF['fusion_id'].tolist()
        flagMultiGeneList=subsetFusDF['flag_multigene'].tolist()
        # Iterate over fusion, starting at the second fusion, and create intron coordinates

        for fus in range(1,len(fusIdList)):
            intronStart=fusStopList[fus-1]+1 #intron starts at the stop position of the previous fusion
            intronStop=fusStartList[fus] #intron stops 1bp before previous fusion
            intronSize=intronStop-intronStart #use this later to delete fusions <1bp in length
            intronChr=chrList[fus]
            intronGene=geneID
            intronFus1=fusIdList[fus-1]
            intronFus2=fusIdList[fus]
            intronID=str(intronFus1) + "|intron|" + str(intronFus2)
            # Keep multigene flags for both introns. Introns made from these fusions will be removed since they
            # are ambiguous
            intronFlag1=flagMultiGeneList[fus-1]
            intronFlag2=flagMultiGeneList[fus]

            #append to intron lists
            intronIDList.append(intronID)
            intronStartList.append(intronStart)
            intronStopList.append(intronStop)
            intronSizeList.append(intronSize)
            intronChrList.append(intronChr)
            intronGeneList.append(intronGene)
            intronFus1List.append(intronFus1)
            intronFus2List.append(intronFus2)
            intronFus1Multigene.append(intronFlag1)
            intronFus2Multigene.append(intronFlag2)

    # Create new dataframe of introns
    intronAllDF = pd.DataFrame(
        {'intron_id': intronIDList, 'chr': intronChrList, 'intron_start': intronStartList,
         'intron_stop': intronStopList, 'intron_size': intronSizeList, 'gene_id': intronGeneList,
         'exonic_region_id_5prime': intronFus1List, 'exonic_region_id_3prime': intronFus2List,
         'flag_multigene_5prime_exon':intronFus1Multigene, 'flag_multigene_3prime_exon':intronFus2Multigene})

    # Drop introns from multigene fusions, as these are ambiguous
    intronAllDF = intronAllDF[intronAllDF['flag_multigene_5prime_exon'] == 0]
    intronAllDF = intronAllDF[intronAllDF['flag_multigene_3prime_exon'] == 0]
    # Drop introns with a 0 bp length
    intronAllDF = intronAllDF[intronAllDF['intron_size'] > 0]

    # Export annotation and BED file for coverage
    intronAllDF = intronAllDF.sort_values(by=['chr', 'intron_start','intron_stop'])
    intronAllDF = intronAllDF[['intron_id', 'chr', 'intron_start','intron_stop','intron_size',
                               'gene_id','exonic_region_id_5prime','exonic_region_id_3prime',
                               'flag_multigene_5prime_exon','flag_multigene_3prime_exon']]
    intronBedDF = intronAllDF[['chr', 'intron_start','intron_stop','intron_id']]
    with open(args.outCSV, 'w') as outFile:
        intronAllDF.to_csv(outFile, encoding='utf-8', index=False)
    with open(args.outBED, 'w') as outFile:
        intronBedDF.to_csv(outFile, encoding='utf-8', sep='\t', header=None, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    # print("Starting script...")
    main()
    logger.info('Script complete')
    # print("Script done!")
