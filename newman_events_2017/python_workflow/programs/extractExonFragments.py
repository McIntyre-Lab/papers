#!/usr/bin/env py27
#######################################################################################################################
#
# DATE: 2018-01-04
# NAME: extractExonFragments.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This program takes a formatted GFF3 file and groups overlapping exons. Each group of exons is then split
# into fragments to identify which fragments ("chunks") are shared by which exons. The output is a CSV that contains
# the exon group start and stop positions (which should be the same coordinates as fusions from buildFusions.py),
# the fragment start and stop positions, and the list of exons assigned to each fragment.
#
# REQUIRED PACKAGES: argparse  (tested with v1.1)
#                    gffutils  (tested with v0.8.5)
#                    logging   (tested with v0.5.1.2)
#
#######################################################################################################################

# Built-in packages
import argparse
import logging

# Add-on packages
import gffutils

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF file")
    parser.add_argument("--output", dest="outCSV", action='store', required=True, help="Output CSV file name")

    args = parser.parse_args()
    return args

def main():
    # Get the database
    logger.info('Connecting to the GFF file')
    db_fn=args.gffInput + '.db'
    db=gffutils.FeatureDB(db_fn)
    
    # Open output file for writing, and write column headers
    csvOut = open(args.outCSV, 'w')
    csvHeader = ','.join(['chr','region_start','region_stop','fragment_start','fragment_stop','exon_id']) + "\n"
    csvOut.write(csvHeader)

    # Get list of chromosomes
    chromosomes=db.features_of_type(['chromosome','chromosome_arm','golden_path_region'])
    
    # Take each chromosome and identify overlapping exons
    for chrom in chromosomes:
        # iterate through exons within chromomsome
        exons = list(db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end), order_by='start'))
        try:
            exons.append(exons[len(exons)-1])
        except:
            continue
        index=1
        exonIndex=0
        lastExon=len(exons)-1 # Check for last exon
        groupExonIDs=[]
        groupExonDesc=[]
        groupExonPos=[]
        groupExons=[]
        groupExonStart=[]
        groupExonStop=[]
        groupStart=[]
        groupStop=[]
        exonGroup=0
        exonNameList=[]
        # First group all overlapping exons into exon groups
        groupExonIDsTemp=[]
        groupExonDescTemp=[]
        groupExonPosTemp=[]
        groupExonsTemp=[]
        groupExonStartTemp=[]
        groupExonStopTemp=[]
        logger.info('Grouping overlapping exons together for chromosome {0}'.format(chrom.id))
        #print("Grouping overlapping exons together for chromosome " + str(chrom.id) + "...")
        for exon in exons:
            # set lists, one for exonIDs, one for feature description (start, end), one for bp position
            if index == 1:
                # Get first exon info -- want start, stop and name
                exonStart=exon.start-1
                exonStop=exon.end
                exonID=''.join(exon.attributes['Name'])
                exonDescStart='start'
                exonDescStop='end'
                currExonStart=exon.start-1
                currExonStop=exon.end
                currExonID=''.join(exon.attributes['Name'])
                currExonDescStart='start'
                currExonDescStop='end'
                index=index+1
                exonIndex=exonIndex+1
                groupExonIDsTemp.append(exonID)
                groupExonIDsTemp.append(exonID)
                groupExonDescTemp.append(exonDescStart)
                groupExonDescTemp.append(exonDescStop)
                groupExonPosTemp.append(exonStart)
                groupExonPosTemp.append(exonStop)
                groupExonsTemp.append(exonID)
                groupExonStartTemp.append(exonStart)
                groupExonStopTemp.append(exonStop)
                groupStart.append(exonStart)
            else:
                # check if current exon overlaps with previous exon, then append to current group
                exonStart=exon.start-1
                exonStop=exon.end
                exonID=''.join(exon.attributes['Name'])
                
                # Overlapping exons check
                # Find pieces of exons that overlap
                if (exonStart <= currExonStop):     # If exon overlaps with current exon
                    # if exon overlaps check to see what has the bigger coordinates
                    if (exonStop >= currExonStop):  # If exon ends later than current exon
                        currExonStart=currExonStart     # Keep current start
                        currExonStop=exonStop           # Make new end point
                        currExonID=exonID
                    exonIndex=exonIndex+1
                    exonDescStart='start'
                    exonDescStop='end'
                    groupExonIDsTemp.append(exonID)
                    groupExonIDsTemp.append(exonID)
                    groupExonDescTemp.append(exonDescStart)
                    groupExonDescTemp.append(exonDescStop)
                    groupExonPosTemp.append(exonStart)
                    groupExonPosTemp.append(exonStop)
                    groupExonsTemp.append(exonID)
                    groupExonStartTemp.append(exonStart)
                    groupExonStopTemp.append(exonStop)
                else:                               # Otherwise start a new group
                    groupExonIDs.append(groupExonIDsTemp)
                    groupExonDesc.append(groupExonDescTemp)
                    groupExonPos.append(groupExonPosTemp)
                    groupExons.append(groupExonsTemp)
                    groupExonStart.append(groupExonStartTemp)
                    groupExonStop.append(groupExonStopTemp)
                    groupStop.append(exonStop)
                    exonGroup=exonGroup+1
                    groupExonIDsTemp=[]
                    groupExonDescTemp=[]
                    groupExonPosTemp=[]
                    groupExonsTemp=[]
                    groupExonStartTemp=[]
                    groupExonStopTemp=[]
                    groupStartTemp=[]
                    groupStopTemp=[]             
                    # Get pieces for next exon
                    exonStart=exon.start-1
                    exonStop=exon.end
                    exonID=''.join(exon.attributes['Name'])
                    exonDescStart='start'
                    exonDescStop='end'
                    currExonStart=exon.start-1
                    currExonStop=exon.end
                    currExonID=''.join(exon.attributes['Name'])
                    currExonDescStart='start'
                    currExonDescStop='end'
                    index=index+1
                    exonIndex=exonIndex+1
                    groupExonIDsTemp.append(exonID)
                    groupExonIDsTemp.append(exonID)
                    groupExonDescTemp.append(exonDescStart)
                    groupExonDescTemp.append(exonDescStop)
                    groupExonPosTemp.append(exonStart)
                    groupExonPosTemp.append(exonStop)
                    groupExonsTemp.append(exonID)
                    groupExonStartTemp.append(exonStart)
                    groupExonStopTemp.append(exonStop)
                    groupStartTemp.append(exonStart)
            if exonIndex == lastExon:
                #Append currently stored temp data to array
                groupExonIDs.append(groupExonIDsTemp)
                groupExonDesc.append(groupExonDescTemp)
                groupExonPos.append(groupExonPosTemp)
                groupExons.append(groupExonsTemp)
                groupExonStart.append(groupExonStartTemp)
                groupExonStop.append(groupExonStopTemp)
                groupStop.append(exonStop)
                exonGroup=exonGroup+1
                groupExonIDsTemp=[]
                groupExonDescTemp=[]
                groupExonPosTemp=[]
                groupExonsTemp=[]
                groupExonStartTemp=[]
                groupExonStopTemp=[]
                groupStartTemp=[]
                groupStopTemp=[]             
                #Get ID start stop parents, etc.
                exonStart=exon.start-1
                exonStop=exon.end
                exonID=''.join(exon.attributes['Name'])
                exonDescStart='start'
                exonDescStop='end'
                currExonStart=exon.start-1
                currExonStop=exon.end
                currExonID=''.join(exon.attributes['Name'])
                currExonDescStart='start'
                currExonDescStop='end'
                index=index+1
                exonIndex=exonIndex+1
        
        # for each group of exons, extra the IDs, feature descriptions and positions and output to new lists
        logger.info('Done grouping overlapping exons for chromosome {0}'.format(chrom.id))
        logger.info('"Dividing exon groups for chromosome {0} into fragments'.format(chrom.id))
        #print("Done grouping overlapping exons for chromosome " + str(chrom.id) +"!")
        #print("Dividing exon groups for chromosome " + str(chrom.id) +" into fragments...")
        for group in range(0,len(groupExonDesc)):
            exonIDs=groupExonIDs[group]
            featureIDs=groupExonDesc[group]
            chromPos=groupExonPos[group]
            # Sort new lists by chromosomal position
            sorted_lists = sorted(zip(exonIDs, featureIDs, chromPos), key=lambda x: x[2])
            exonIDs, featureIDs, chromPos = [[x[i] for x in sorted_lists] for i in range(3)]
            # iterate through features within exon group
            itemIndex=0
            fragIndex=0
            fragExonList=[]
            fragFeatureList=[]
            fragStartList=[]
            fragStopList=[]
            fragExons=[]
            fragFeatures=[]
            # Try appending last fragment again, so that when I hit the last piece I can output whatever is currently in memory
            #featureIDs.append(featureIDs[len(featureIDs)-1])
            # set the count for the last feature
            lastFeature=len(featureIDs)
            for feat in range(0,len(featureIDs)):
                if itemIndex == 0:
                    # get first info in lists
                    currPos=chromPos[feat]
                    currFeat=featureIDs[feat]
                    currExon=exonIDs[feat]
                    prevPos=chromPos[feat]
                    prevFeat=featureIDs[feat]
                    prevExon=exonIDs[feat]
                    fragExons.append(currExon)
                    fragFeatures.append(currFeat)
                    fragStart=prevPos
                    itemIndex=itemIndex+1
                else:
                    try:
                        currPos=chromPos[feat]
                        currFeat=featureIDs[feat]
                        currExon=exonIDs[feat]
                    except:
                        continue
                    if currPos == prevPos:
                        fragExons.append(currExon)
                        fragFeatures.append(currFeat)
                        continue
                    else:
                        fragStop=currPos
                        fragExons.append(currExon)
                        fragFeatures.append(currFeat)
                        fragExonList.append(fragExons)
                        fragFeatureList.append(fragFeatures)
                        fragStartList.append(fragStart)
                        fragStopList.append(fragStop)
                        itemIndex=itemIndex+1
                        # Write current chunk to output
                        fragExonOut='|'.join(fragExons)
                        fragFeaturesOut='|'.join(fragFeatures)
                        fragNum=fragIndex+1
                        # Start new fragment
                        fragExons=[]
                        fragFeatures=[]
                        try:
                            prevPos=chromPos[feat]
                            prevFeat=featureIDs[feat]
                            prevExon=exonIDs[feat]
                        except:
                            continue
                        fragStart=currPos
                        fragIndex=fragIndex+1
                if itemIndex == lastFeature:
                    fragStop=currPos
                    fragExons.append(currExon)
                    fragFeatures.append(currFeat)
                    fragExonList.append(fragExons)
                    fragFeatureList.append(fragFeatures)
                    fragStartList.append(fragStart)
                    fragStopList.append(fragStop)
                    #Write current lists to output
                    fragExonOut='|'.join(fragExons)
                    fragFeaturesOut='|'.join(fragFeatures)
                    fragNum=fragIndex+1
            # iterate through exons, append exonIDs to list
            fragExonIDList=[]
            for cx in range(0,len(fragStartList)):
                fragTotalStart=fragStartList[0]
                fragTotalStop=fragStopList[len(fragStopList)-1]
                fragExonIDList2=[]
                for ex in range(0,len(groupExons)):
                    for ey in range(0,len(groupExons[ex])):
                        ind=1
                        if groupExonStart[ex][ey] <= fragStartList[cx] and groupExonStop[ex][ey] >= fragStopList[cx]:
                            if ind == 1 :
                                fragExonIDList2.append(groupExons[ex][ey])
                                ind=ind+1
                            elif fragExonIDList2[len(fragExonIDList2)-1] != groupExons[ex][ey]:
                                fragExonIDList2.append(groupExons[ex][ey])
                if fragStartList[cx] != fragStopList[cx]:
                    # Avoid adding in fragments of size 0
                    csvEntry = ','.join([str(chrom.id),str(fragTotalStart),str(fragTotalStop),str(fragStartList[cx]),str(fragStopList[cx]),str('|'.join(fragExonIDList2))]) + "\n"
                    csvOut.write(csvEntry)

        logger.info('Done dividing exon groups for chromosome {0} into fragments'.format(chrom.id))
        #print("Done dividing exon groups for chromosome " + str(chrom.id) + " into fragments!")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    #print("Starting script...")
    main()
    logger.info('Script complete')
    #print("Script done!")
