#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-04-02
# NAME: extract_junction_sequence.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script takes the coordinate-collapsed junction annotations and extracts FASTA sequences for each.
# It then outputs four files: junction annotations, unique junction sequences, sequence-to-junction index, BED file
# for calculating junction read coverage
#
# REQUIRED PACKAGES: argparse   (tested with v1.1)
#                    pandas     (tested with v0.19.2)
#                    logging    (tested with v0.5.1.2)
#                    pybedtools (tested with v0.7.10)
#
#######################################################################################################################


# Packages
import argparse
import logging

import pandas as pd
import pybedtools
import sys
from io import StringIO



def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import fragments and fusion annotation files")
    parser.add_argument("--input-junction-file", dest="inJuncInfo", action='store', required=True,
                        help="Input junction annotation CSV file with donor/acceptor information, collapsed by"
                             "genome coordinate")
    parser.add_argument("--input-fasta-file", dest="inFasta", action='store', required=True,
                        help="Input genome FASTA file, containing all chromosome/contig sequences")
    parser.add_argument("--output-junction-annotation", dest="outJuncAnnot", action='store', required=True,
                        help="Output junction annotation CSV file with donor/acceptor information, collapsed by"
                             "genome coordinate")
    parser.add_argument("--output-junction-to-seq-index", dest="outSeqIndex", action='store', required=True,
                        help="Output annotation file that links junctions to unique sequences")
    parser.add_argument("--output-junction-sequences", dest="outJuncSeq", action='store', required=True,
                        help="Output multi-sequence FASTA file containing all unique junction sequences")
    parser.add_argument("--output-coverage-bed", dest="outJuncBED", action='store', required=True,
                        help="Output 3-column BED-file for coverage counts")

    args = parser.parse_args()
    return args

def main():
    # Import collapsed junction annotations
    juncInDF = pd.read_csv(args.inJuncInfo, sep=",", 
                           usecols=['junction_id','chr','donor_start','donor_stop','acceptor_start','acceptor_stop',
                                  'strand','transcript_id','gene_id','annotation_frequency', 'flag_multigene',
                                    'flag_junction_annotated','flag_border_junction', 'flag_alt_donor',
                                    'flag_alt_acceptor', 'flag_exonskip'])

    # Set genome FASTA file reference
    fasta_in=args.inFasta
    # Set up lists for appending
    juncIDlist=[]
    seqList=[]
    bedList=[]
    for indexJunc, rowJunc in juncInDF.iterrows():
        idJunc=str(rowJunc['junction_id'])
        chrJunc=rowJunc['chr']
        donorStartJunc=rowJunc['donor_start']
        donorStopJunc=rowJunc['donor_stop']
        acceptorStartJunc=rowJunc['acceptor_start']
        acceptorStopJunc=rowJunc['acceptor_stop']
        strandJunc=str(rowJunc['strand'])
        flagBorderJunc=rowJunc['flag_border_junction']
        # For each junction, create a BED-12 entry as a space-delimited string, as junctions consist of a donor and
        # an acceptor sequence
        # BED12: chr, start, stop, ID, score, strand, start, stop, rbg, num_blocks, block_size, block_start
        # Border junctions can be considered a single block of sequence, so these need to treated slightly different
        if flagBorderJunc == 1:
            blockCount= "1"
            blockStart = "0"
            blockSize = str((acceptorStopJunc - donorStartJunc) + 1)
            if strandJunc == "+":
                totalStart=str(donorStartJunc)
                totalStop=str(acceptorStopJunc+1)
            else:
                totalStart = str(donorStartJunc-1)
                totalStop = str(acceptorStopJunc)
        else:
            blockCount="2"
            totalStart=str(donorStartJunc)
            totalStop = str(acceptorStopJunc)
            blockSize=str(donorStopJunc-donorStartJunc)+","+str(acceptorStopJunc-acceptorStartJunc)
            blockStart="0," + str(acceptorStartJunc-donorStartJunc)
        bedEntryStr=str(chrJunc) + " " + totalStart + " " + totalStop + " " + idJunc + " . " + strandJunc + " " + \
                    totalStart + " " + totalStop + " 0,0,0 " + blockCount + " " + blockSize + " " + blockStart
        bedList.append(bedEntryStr)
    # Parse string into pybedtools.BedTool to set up coordinate and extract sequences
    bedListStr='\n'.join(bedList)
    bedEntry = pybedtools.BedTool(bedListStr, from_string=True)
    # Extract sequence, allow for strandness and splitting
    bedSeq = bedEntry.sequence(fi=fasta_in, s=True, name=True, split=True, tab=True)
    # Get junction sequence
    #juncSeq = str(open(bedSeq.seqfn).read()).split("\n")[1].upper()
    juncSeq = StringIO(str(open(bedSeq.seqfn).read()))
    # Append to lists
    junc2seqDF = pd.read_table(juncSeq, sep="\t", header=None, names=['junction_id','sequence'])
    junc2seqDF['junction_id'] = junc2seqDF['junction_id'].str.replace(r"\(\-\)|\(\+\)", "")
    uniqSeqList=junc2seqDF['sequence'].sort_values().drop_duplicates(keep='first').tolist()
    seqLenList =[]
    uniqIdList =[]
    uniqNum = 1
    for seq in range(0,len(uniqSeqList)):
        seqLen=len(uniqSeqList[seq])
        uniqSeqID="junction_" + str(uniqNum)
        uniqNum=uniqNum+1
        seqLenList.append(seqLen)
        uniqIdList.append(uniqSeqID)
    uniqSeqDF = pd.DataFrame({'sequence_id': uniqIdList, 'sequence':uniqSeqList, 'sequence_length': seqLenList})
    # Merge unique sequences to junctions and create a junction-to-sequence index
    uniqSeqDF = uniqSeqDF.sort_values(by=['sequence'])
    junc2seqDF = junc2seqDF.sort_values(by=['sequence'])
    junc2uniqDF=pd.merge(junc2seqDF, uniqSeqDF, how='inner', on=['sequence'])
    junc2uniqDF=junc2uniqDF[['junction_id','sequence_id','sequence']]
    # Create 3-column BED file for coverage: seqID, start(0), stop(length)
    bedCoverageDF=uniqSeqDF[['sequence_id','sequence_length']]
    bedCoverageDF['start'] = pd.Series(uniqSeqDF['sequence_length']-uniqSeqDF['sequence_length'])
    bedCoverageDF = bedCoverageDF[['sequence_id','start','sequence_length']]
    # Write outputs
    ## BED file
    with open(args.outJuncBED, 'w') as outFile:
        bedCoverageDF.to_csv(outFile, encoding='utf-8', index=False, sep="\t", header=False)
    ## Junction to sequence index
    with open(args.outSeqIndex, 'w') as outFile:
        junc2uniqDF.to_csv(outFile, encoding='utf-8', index=False)
    ## Annotation CSV
    juncAnnotDF = juncInDF[['junction_id','chr','donor_start','donor_stop','acceptor_start','acceptor_stop','strand',
                            'transcript_id','gene_id','annotation_frequency','flag_multigene','flag_junction_annotated',
                            'flag_border_junction','flag_alt_donor','flag_alt_acceptor','flag_exonskip']]
    with open(args.outJuncAnnot, 'w') as outFile:
        juncAnnotDF.to_csv(outFile, encoding='utf-8', index=False)
    ## FASTA file of sequences
    fastaOut=[]
    for indexSeq, rowSeq in uniqSeqDF.iterrows():
        seqName=">"+str(rowSeq['sequence_id'])
        seqSeq=rowSeq['sequence']
        fastaOut.append(seqName)
        fastaOut.append(seqSeq)
    with open(args.outJuncSeq, 'w') as outputFile:
        outputFile.write("\n".join(str(i) for i in fastaOut))

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    main()
    logger.info('Script complete')
