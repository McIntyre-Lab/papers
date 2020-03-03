#!/usr/bin/env python
#Standart Libraries
from itertools import groupby
import collections
import argparse 
import logging
import re
import os

#AddOn Libraries
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#Args options
def getOptions():
    """Function to pull in arguments"""
    parser = argparse.ArgumentParser(description="Takes a single-end (SE) file and splits out unique and duplicate reads.")
    parser.add_argument("-r1", dest="r1", action='store', required=True, help="Name of read 1 or SE FASTQ file [required]")
    parser.add_argument("-r2", dest="r2", action='store', required=False, help="Name of read 2 or PE FASTQ file [required]",default=False)
    parser.add_argument("--outdir", dest="odir", action='store', required=True, help="Directory to store FASTQ output files [required]")
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Path and name of the sumamary table of read counts in csv format [require]")
    parser.add_argument("-t", "--table", dest="tname", action='store', required=False, default=False, help="Path and name of the table showing each sequence and the number of reads with that sequence in tsv format [Optional]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file")
    args =parser.parse_args()
    return (args)
    
def readFastQ_SE(fastq_path):
    fastq_file = open(fastq_path,'r')
    fastq_generator=FastqGeneralIterator(fastq_file)
    readDict = {re.sub('/[1-2]','',header).split(' ')[0]:(seq,qual) for header,seq,qual in fastq_generator}
    readTupl = [(header,readDict[header][0]) for header in readDict]
    return (readDict,readTupl)

def readFastQ_PE(fastq_path1,fastq_path2):
    fastq_file1 = open(fastq_path1,'r')
    fastq_file2 = open(fastq_path2,'r')

    fastq_generator1=FastqGeneralIterator(fastq_file1)
    fastq_generator2=FastqGeneralIterator(fastq_file2)

    readDict1 = {re.sub('/[1-2]','',header).split(' ')[0]:(seq,qual) for header,seq,qual in fastq_generator1}
    readDict2 = {re.sub('/[1-2]','',header).split(' ')[0]:(seq,qual) for header,seq,qual in fastq_generator2}

    readTupl = []
    #in 1 but not in 2
    readUnprd12 = []
    #in 2 but not in 1
    readUnprd21 = []

    for header in readDict1:
    	if header in readDict2:
    		readTupl.append((header,readDict1[header][0]+readDict2[header][0]))
    	else:
    		readUnprd12.append((header,readDict1[header][0]))

    for header in readDict2:
    	if header in readDict1: continue
    	else:
    		readUnprd21.append((header,readDict2[header][0]))

    if len(readUnprd12):
    	readUnprd12=False
    if len(readUnprd21):
    	readUnprd21=False

    return (readDict1,readDict2,readTupl,readUnprd12,readUnprd21)

def getHeadersList(readTupl):
	disHead=[]
	dupHead=[]
	sduHead=[]
	uniHead=[]

	flag_inDis_off=True
	for i in range(len(readTupl)):
		if i!=0 and i<(len(readTupl)-1):
			if readTupl[i][1]==readTupl[i+1][1]:
				dupHead.append(readTupl[i][0])
				if flag_inDis_off:
					disHead.append(readTupl[i][0])
					sduHead.append(readTupl[i][0])
					flag_inDis_off=False
			elif readTupl[i][1]==readTupl[i-1][1]:
				dupHead.append(readTupl[i][0])
				flag_inDis_off=True
			else:
				disHead.append(readTupl[i][0])
				
		elif i!=0 and i+1==len(readTupl):
			if readTupl[i][0]==readTupl[i-1][0]:
				dupHead.append(readTupl[i][0])
			else:
				disHead.append(readTupl[i][0])
				
		else:
			disHead.append(readTupl[i][0])
			if readTupl[i][0]==readTupl[i+1][0]:
				dupHead.append(readTupl[i][0])
				sduHead.append(readTupl[i][0])
				flag_inDis_off=False
	
	uniHead=list(set(disHead)-set(sduHead))
	return (disHead,dupHead,sduHead,uniHead)

#Print ut CVS format file.
def writeTableSE (UNI,DIS,DUP,SDU,oname):
	total = len(DIS)+len(DUP)-len(SDU)
	
	percent_uniq_seq = float(len(UNI))/float(total) * 100
	percent_uniq_reads = float(len(UNI))/float(total) * 100
	percent_dist_reads = float(len(DIS))/float(total) * 100
	percent_dup_reads = float(len(DUP))/float(total) * 100
	
	ON = open(oname,'w')
	ON.write("file_name,total_reads,num_unique_reads,per_unique_reads,num_distinct_reads,per_distinct_reads,num_duplicate_reads,per_duplicate_reads\n")
	out=[os.path.basename(oname),total,len(UNI),percent_uniq_reads,len(DIS),percent_dist_reads,len(DUP),percent_dup_reads]
	ON.write(','.join([str(x) for x in out])+"\n")
	ON.close()

def writeOutput (headList,readDict,readNumber,OUTFILE):
	for head in headList:
		OUTFILE.write ('\n'.join(['@'+head+readNumber,readDict[head][0],'+',readDict[head][1],'']))
		
def readsCount (readTupl,tname):
	tname = os.path.abspath(tname)
	
	justReads = [read for ID,read in readTupl]
	seqCounter = collections.Counter(justReads)
	seqCounter = [(seqCounter[seq],seq)for seq in seqCounter]
	seqCounter.sort(key=lambda x: x[0], reverse=True)

	TFILE = open(tname,'w')
	for seq in seqCounter:
		print >> TFILE,'\t'.join([str(seq[0]),seq[1]])
	TFILE.close()

def main(args):
	fastq_path1 = args.r1
	fastq_path2 = args.r2
	
	if args.r2:
		if args.log:
			logger.info("Reading '%s' '%s' " % (args.r1, args.r2))
		readDict1,readDict2,readTupl,readUprd12,readUprd21 = readFastQ_PE(fastq_path1,fastq_path2)
		if args.log:
			logger.info("Finished reading '%s' '%s' " % (args.r1, args.r2))
		readTupl.sort(key=lambda x: x[1])
		readUprd12.sort(key=lambda x: x[1])
		readUprd21.sort(key=lambda x: x[1])

		disHead,dupHead,sduHead,uniHead = getHeadersList(readTupl)

		odir = os.path.abspath(args.odir)

		bname1 = os.path.basename(args.r1)
		name1 = os.path.splitext(bname1)[0]
		bname2 = os.path.basename(args.r2)
		name2 = os.path.splitext(bname2)[0]

		DIS1 = open (os.path.join(odir, name1 + '_distinct.fq'),'w')
		DUP1 = open (os.path.join(odir, name1 + '_duplicate.fq'),'w')
		SDU1 = open (os.path.join(odir, name1 + '_single_duplicate.fq'),'w')
		UNI1 = open (os.path.join(odir, name1 + '_uniq.fq'),'w')
		DIS2 = open (os.path.join(odir, name2 + '_distinct.fq'),'w')
		DUP2 = open (os.path.join(odir, name2 + '_duplicate.fq'),'w')
		SDU2 = open (os.path.join(odir, name2 + '_single_duplicate.fq'),'w')
		UNI2 = open (os.path.join(odir, name2 + '_uniq.fq'),'w')

		writeOutput(disHead, readDict1, '/1', DIS1)
		writeOutput(dupHead, readDict1, '/1', DUP1)
		writeOutput(sduHead, readDict1, '/1', SDU1)
		writeOutput(uniHead, readDict1, '/1', UNI1)
		writeOutput(disHead, readDict2, '/2', DIS2)
		writeOutput(dupHead, readDict2, '/2', DUP2)
		writeOutput(sduHead, readDict2, '/2', SDU2)
		writeOutput(uniHead, readDict2, '/2', UNI2)
		writeTableSE(uniHead,disHead,dupHead,sduHead,args.oname)

		if readUprd12:
			DISU1 = open (os.path.join(odir, name1 + '_unpaired_distinct.fq'),'w')
			DUPU1 = open (os.path.join(odir, name1 + '_unpaired_duplicate.fq'),'w')
			SDUU1 = open (os.path.join(odir, name1 + '_unpaired_single_duplicate.fq'),'w')
			UNIU1 = open (os.path.join(odir, name1 + '_unpaired_uniq.fq'),'w')
			writeOutput(disHead, readDict1, '/1', DISU1)
			writeOutput(dupHead, readDict1, '/1', DUPU1)
			writeOutput(sduHead, readDict1, '/1', SDUU1)
			writeOutput(uniHead, readDict1, '/1', UNIU1)
		if readUprd21:
			DISU2 = open (os.path.join(odir, name2 + '_unpaired_distinct.fq'),'w')
			DUPU2 = open (os.path.join(odir, name2 + '_unpaired_duplicate.fq'),'w')
			SDUU2 = open (os.path.join(odir, name2 + '_unpaired_single_duplicate.fq'),'w')
			UNIU2 = open (os.path.join(odir, name2 + '_unpaired_uniq.fq'),'w')
			writeOutput(disHead, readDict2, '/2', DISU2)
			writeOutput(dupHead, readDict2, '/2', DUPU2)
			writeOutput(sduHead, readDict2, '/2', SDUU2)
			writeOutput(uniHead, readDict2, '/2', UNIU2)
		if args.tname:
			tname=args.tname
			readsCount(readTupl,tname)


	else:
		if args.log:
			logger.info("Reading '%s' " % (args.r1))
		readDict,readTupl = readFastQ_SE(fastq_path1)
		if args.log:
			logger.info("Finished reading '%s' " % (args.r1))
		readTupl.sort(key=lambda x: x[1])
		disHead,dupHead,sduHead,uniHead = getHeadersList(readTupl)

		odir = os.path.abspath(args.odir)
		bname = os.path.basename(args.r1)
		name = os.path.splitext(bname)[0]

		DIS = open (os.path.join(odir, name + '_distinct.fq'),'w')
		DUP = open (os.path.join(odir, name + '_duplicate.fq'),'w')
		SDU = open (os.path.join(odir, name + '_single_duplicate.fq'),'w')
		UNI = open (os.path.join(odir, name + '_uniq.fq'),'w')

		writeTableSE(uniHead,disHead,dupHead,sduHead,args.oname)
		writeOutput(disHead, readDict, '', DIS)
		writeOutput(dupHead, readDict, '', DUP)
		writeOutput(sduHead, readDict, '', SDU)
		writeOutput(uniHead, readDict, '', UNI)
	
if __name__=='__main__':
	args = getOptions()
	if args.log:
		logging.basicConfig(filename=(os.path.abspath(args.log)),level=logging.DEBUG)
		logger = logging.getLogger()
	main(args)
	if args.log:
		logger.info("Script complete.")
