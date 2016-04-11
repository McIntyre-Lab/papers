#!/usr/bin/env python
#########################################################################################################
# SCRIPT: windowMaker_mai1.py
# 
# AUTHOR: Miguel Ibarra Arellano (miguelib@ufl.edu)
# 
# DESCRIPTION: This script takes a Variant Call File (VCF) and splits it into windows of a user-defined
# size. Two inputs are required: the VCF file of variants you want to create windows for, and a design
# file that defines window sizes (tab separated or comma separated. This file should have the headings:
#     window_name      (what the windows will be labelled)
#     window_size      (size of windows in kilobases (kb))
#     chrom            (chromosome)
#     windowing_start  (position to start windowing)
#     windowing_end    (position to stop windowing)
# 
# The output is a series of VCF files for each window, and a summary file with a list of all the Windows
# and the distribution of variants per file.
#
#########################################################################################################

#Standard Libraries
import argparse
import os
import logging
import itertools
import shutil

#AddOn Libraries
import numpy as np
import multiprocessing as mp

#Getting all the arguments
def getOptions():
	"""Function to pull arguments"""
	parser = argparse.ArgumentParser(description="Takes a file with m/z and RT and compares per feature against another file with m/z and RT")

	required = parser.add_argument_group(title='Requiered Input', description='Requiered input to the program')
	required.add_argument('-v','--VCF',dest="VCF", action="store",required=True,help="Name of the VCF file")
	required.add_argument('-e','--especFILE',dest="especFILE", action="store",required=True,help="Name of the first file (could be either csv or tsv format)")
	required.add_argument('-o','--OUTDIR',dest="OUTDIR",action="store",required=True,help="Name of the Output Directory")

	optional = parser.add_argument_group(title="Optional Input", description="Optional Input to the program")
	optional.add_argument('-g',"--LOG",dest="LOG",action="store",required=False,help="Path and name of the log file")

	args = parser.parse_args()
	return(args);

#readTable, reads any table foramated file from a path and returns a numpy array  with it.
#Input:
#	filePath:	path of the file
#Output:
#	tabledFile:	numpy array objet with the file. [[],[],[],[],[]]
def readTable (filePath):
	"""Function to read Files"""
	bname = os.path.basename(filePath)
	name,format = os.path.splitext(filePath)

	if args.LOG:logger.info("Reading File '%s' " %(filePath))
	if format == '.csv':
		tabledFile = np.genfromtxt(fname=filePath, skip_header=1, dtype='str', delimiter=',',comments='#')
	else:
		tabledFile = np.genfromtxt(fname=filePath, skip_header=1, dtype='str', delimiter='\t', comments='#')
	if args.LOG:logger.info("Finished reading File '%s' " %(filePath))
	return (tabledFile)

#fixedWindowMaker, creates windows from a given VCF file following the especifications in especifications file.
#Input:
#	especFILE:	numpy array objet with the file. [[],[],[],[],[]]
#	VCFFile		path of the VCF file
#	ODIR:		path to the output directory
#Output:
#	tabledFile:	numpy array objet with the file. [[],[],[],[],[]]
def fixedWindowMaker (especFILE,VCFFile,ODIR):
	#Solves litle bug when the especFILE 
	if str(type(especFILE[0]))=='<type \'numpy.string_\'>':
		especFILE=[especFILE]
 
	WINDOWLIST = open(os.path.join(ODIR,"widowList.tsv"),'w')
	print >> WINDOWLIST,"WindowName\tWindowSize\tChrom\tWindowStart\tWindowEnd"

	for spec in especFILE:
		windowingEnd = int(spec[4])
		windowStart = int(spec[3])
		windowSize = int(float(spec[1])*1000)
		windowName = spec[0]
		chromName = spec[2]

		windowEnd = windowSize
		count = 0

		WINDOwDIR = os.path.join(ODIR,"chrom_"+chromName)
		if os.path.exists(WINDOwDIR):
			shutil.rmtree(WINDOwDIR)
		os.makedirs(WINDOwDIR)
		
		while windowStart < windowingEnd:
			windowFileName = windowName+'_'+spec[1]+"kb_chrom"+chromName+"_"+str(count)
			WINDOWFILE = open(os.path.join(WINDOwDIR, windowFileName+'.vcf'),'w')
			
			VCFFile = open(args.VCF,'r')#Opens the file and goes through line by line
			numVariants = 0
			for VCFLine in VCFFile:
				if VCFLine.startswith("##"):continue #Ignores the comentaries
				VCFLine = VCFLine.strip()
				if VCFLine.startswith("#"): #Prints the header
					print >> WINDOWFILE,VCFLine
					continue
				VCFLine = VCFLine.split('\t') #Splits line by tabs
				if ((int(VCFLine[1]) >= windowStart) and (int(VCFLine[1]) < windowEnd)):#If the variant is into the window then saves the variant in its respective window file
					print >> WINDOWFILE,"\t".join(VCFLine)
					numVariants += 1
			VCFFile.close()
			
			#Print information of the window in the WINDOWLIST file
			print >> WINDOWLIST,'\t'.join(map(str,[windowFileName+".vcf",numVariants,chromName,windowStart,windowEnd]))

			#Calculates new window 
			windowStart = windowEnd+1
			windowEnd = windowEnd+windowSize
			count +=1

def main(args):
	"""Function to call all other functions"""
	especFILE = readTable(args.especFILE)
	ODIR = args.OUTDIR
	fixedWindowMaker(especFILE,args.VCF,ODIR)

if __name__=='__main__':
	args = getOptions()
	if args.LOG:
		logging.basicConfig(filename=(os.path.abspath(args.LOG)),level=logging.DEBUG)
		logger = logging.getLogger()
	main(args)
	if args.LOG:logger.info("Script complete.")
