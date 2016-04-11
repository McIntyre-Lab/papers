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
	required.add_argument('-o','--OUTDIR',dest="OUTDIR",action="store",required=True,help="Name of the Output Directory")

	method = parser.add_mutually_exclusive_group(required=True)
	method.add_argument("-e","--specFILE", dest="specFILE", action="store", help="Path to the spec FILE.")
	method.add_argument("-b","--bedFILE", dest="bedFILE", action="store", help="Path to the .bed File")

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

#fixedWindowMaker, creates windows from a given VCF file following the specifications in specifications file.
#Input:
#	specFILE:	numpy array objet with the file. [[],[],[],[],[]]
#	VCFFile		path of the VCF file
#	ODIR:		path to the output directory
#Output:
#	tabledFile:	numpy array objet with the file. [[],[],[],[],[]]
def fixedWindowMaker (specFILE,VCFFile,ODIR):
	#Solves litle bug when the specFILE 
	if str(type(specFILE[0]))=='<type \'numpy.string_\'>':
		specFILE=[specFILE]
 
	WINDOWLIST = open(os.path.join(ODIR,"widowList.tsv"),'w')
	print >> WINDOWLIST,"WindowName\tWindowSize\tChrom\tWindowStart\tWindowEnd"

	for spec in specFILE:
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
				if ((int(VCFLine[1]) >= windowStart) and (int(VCFLine[1]) < windowEnd)):#If the variant is into the window then saves the variant in its rspective window file
					print >> WINDOWFILE,"\t".join(VCFLine)
					numVariants += 1
			VCFFile.close()
			
			#Print information of the window in the WINDOWLIST file
			print >> WINDOWLIST,'\t'.join(map(str,[windowFileName+".vcf",numVariants,chromName,windowStart,windowEnd]))

			#Calculates new window 
			windowStart = windowEnd+1
			windowEnd = windowEnd+windowSize
			count +=1


def bedWindowMaker (bedTable,VCFFile,ODIR):
	"""
	#bedWindowMaker, creates windows from a given VCF file following the specifications in a bed file.
	#Input:
	#	bedTable:	numpy array objet with the file. [[],[],[],[],[]]
	#	VCFFile		path of the VCF file
	#	ODIR:		output path
	#Output:
	#	tabledFile:	numpy array objet with the file. [[],[],[],[],[]]
	"""

	#Solves litle bug when the specFILE 
	if str(type(bedTable[0]))=='<type \'numpy.string_\'>':
		bedTable=[bedTable]
 
	#Open a WINDOWLIST file, it contains all the necesary information per file
	WINDOWLIST = open(os.path.join(ODIR,"windowList.tsv"),'a')

	#Print header
	print >> WINDOWLIST,"WindowName\tVariants\tChrom\tWindowStart\tWindowEnd"

	#Read VCF file
	VCFFile = np.genfromtxt(fname=args.VCF, skip_header=0,dtype='str',delimiter='\t',comments='##')

	#Creates new file for line in bedTable.
	for bedLine in bedTable:

		#Takes some arguments of bed format https://genome.ucsc.edu/FAQ/FAQformat.html
		chrom = bedLine[0]
		chromStart = int(bedLine[1])
		chromEnd = int(bedLine[2])
		name = bedLine[3]
								
		#Creates an output path for each chromosome
		WINDOwDIR = os.path.join(ODIR,"chrom_"+chrom)
		out=[]

		#Iterates over opened VCF file
		for VCFLine in VCFFile:
			#Saves header,
			if VCFLine[0].startswith("#"):
				out.append("\t".join(VCFLine))
				continue

			###NOTE###
			#This piece of code asumes that VCF file is sorted, if is not then the script may not work

			#Saves variant (gen, chrom specific)
			if (int(VCFLine[1])<=chromStart):
				 continue
			elif((int(VCFLine[1])>=chromStart) and (int(VCFLine[1])<=chromEnd) and (VCFLine[0]==chrom)):
				out.append("\t".join(VCFLine))
				#Creates new folder
				try:
					os.makedirs(WINDOwDIR)
				except:
					pass
			else:
				break

		#Verifies if any variant found
		if len(out)>1:
			print "entre"

			#Opens output VCF file
			WINDOWFILE = open(os.path.join(WINDOwDIR, name+'.vcf'),'w')

			#Write output
			print >> WINDOWFILE,"\n".join(out)

			#Close output VCF file
			WINDOWFILE.close()

		#Print information of the window in the WINDOWLIST file
		print >> WINDOWLIST,'\t'.join(map(str,[chrom+"_"+name+".vcf",len(out)-1,chrom,chromStart,chromEnd]))
	WINDOWLIST.close


def main(args):
	"""Function to call all other functions"""

	#Reading specFILE or BED file
	if (args.bedFILE):
		bedFILE = readTable(args.bedFILE)
	else:
		specFILE = readTable(args.espcFILE)


	#Windowing for spectFILE
	if (args.bedFILE):
		bedWindowMaker(bedFILE,args.VCF,args.OUTDIR)
	else:
		fixedWindowMaker(specFILE,args.VCF,args.OUTDIR)


if __name__=='__main__':
	args = getOptions()
	if args.LOG:
		logging.basicConfig(filename=(os.path.abspath(args.LOG)),level=logging.DEBUG)
		logger = logging.getLogger()
	main(args)
	if args.LOG:logger.info("Script complete.")
