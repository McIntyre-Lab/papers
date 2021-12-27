#!/usr/bin/env python

# Take in BBmerge merged fastq output and split reads into >=50% original read length and <50%
# if both read pairs are empty - delete
# if only 1 read pair is empty - output 

from Bio import SeqIO
import argparse
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Take in BBmerge merged fastq output and split reads into >=(50% original read length + 14bp) and <(50% + 14bp)")

    # Input arguments
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Merged FASTQ file from BBmerge")
    parser.add_argument("-r1", "--R1", dest="R1File", required=True, help="R1 FASTQ file from BBmerge")
    parser.add_argument("-r2", "--R2", dest="R2File", required=True, help="R2 FASTQ file from BBmerge")
    parser.add_argument("-rl", "--readLength", dest="readLength", required=True, type=int, help="Original read length")

    # Output arguments
    parser.add_argument("-g", "--greaterEqualDir", dest="GEdir", required=True, help="Path to directory for BBmerged reads that are greater than or equal to 50% the original read length + 14bp")
    parser.add_argument("-l", "--lessThanDir", dest="Ldir", required=True, help="Path to directory for BBmerged reads that are less than 50% the original read length + 14bp")    
    parser.add_argument("-o", "--widowDir", dest="widowDir", required=True, help="Path to directory for all unmerged R1 and R2 reads that are >= 50% of the readLength + 14bp. Note - Widow reads are reads (R1 or R2) that had and empty/short pair.")

    args = parser.parse_args()
    return args


def main():
    # Check that input file and output directories exist
    if not os.path.isfile(args.inFile) or not os.path.isfile(args.R1File) or not os.path.isfile(args.R2File) or not os.path.isdir(args.GEdir) or not os.path.isdir(args.Ldir) or not os.path.isdir(args.widowDir):
        raise FileNotFoundError

    # Get filename for output files of the same name in the new output directories
    mergeName = os.path.basename(args.inFile)

    R1name =  os.path.basename(os.path.splitext(args.R1File)[0])
    R2name =  os.path.basename(os.path.splitext(args.R2File)[0])
    
    # Set cutoff length
    cutoff = ((args.readLength/2) + 14)
    print("Merged Read cutoff (50% RL) = " + str(cutoff))

    # Import fastq input file and split by sequence length
    listGE = []
    listL = []
    for seqRecord in SeqIO.parse(args.inFile, "fastq"):
        if len(seqRecord) >= cutoff:
            listGE.append(seqRecord)
        else :
            listL.append(seqRecord)
    print("Count of bbmerged sequences found >= 50% RL + 14bp = " + str(len(listGE)))
    print("Count of bbmerged sequences found < 50% RL + 14bp = " + str(len(listL)))

    SeqIO.write(listGE, args.GEdir + "/" + mergeName, "fastq")
    SeqIO.write(listL, args.Ldir + "/" + mergeName, "fastq")

    listR1empty = []
    listR2empty = []
    listR1short = []
    listR2short = []
    listBoth = []
        
    ## for R1 files
    for seqRecord in SeqIO.parse(args.R1File,"fastq"):
        if len(seqRecord) == 0:
            listR1empty.append(seqRecord.id)
        if len(seqRecord) > 0 and len(seqRecord) < cutoff:
            listR1short.append(seqRecord.id)
    print("number R1 reads empty = " + str(len(listR1empty)))
    print("number R1 reads short = " + str(len(listR1short)))

    ## for R2 files
    for seqRecord in SeqIO.parse(args.R2File,"fastq"):
        if len(seqRecord) == 0:
            listR2empty.append(seqRecord.id)    
        if  len(seqRecord) > 0 and len(seqRecord) < cutoff:
            listR2short.append(seqRecord.id)
    print("number R2 reads empty = " + str(len(listR2empty)))
    print("number R2 reads short = " + str(len(listR2short)))

    ## combine lists 
    listComboR1 = listR1empty + listR1short
    listComboR2 = listR2empty + listR2short

    print("# R1 empty + short = " + str(len(listComboR1)))
    print("# R2 empty + short = " + str(len(listComboR2)))

    ## print empty and short readIDs to text files
    f = open(args.widowDir + "/" + R1name + "_list_short_empty.txt", "w+")
    f.writelines("%s\n" % read for read in listComboR1)
    f = open(args.widowDir + "/" + R2name + "_list_short_empty.txt", "w+")
    f.writelines("%s\n" % read for read in listComboR2)


    # find intersection between lists - these are read pairs that are empty or short in both - remove 
    listBoth = list(set(listComboR1) & set(listComboR2))
    print("empty or short read PAIRS = " + str(len(listBoth))) 

            
    ## if R1 read is short or empty then remove from both R1 and R2 paired files
    ## if R2 read is short or empty then remove from both R1 and R2 paired files  
    R1_dict = SeqIO.to_dict(SeqIO.parse(args.R1File, "fastq"))
    R2_dict = SeqIO.to_dict(SeqIO.parse(args.R2File, "fastq"))
    print("total number in R1 = " + str(len(R1_dict)))
    print("total number in R2 = " + str(len(R2_dict)))
    for key in listBoth:
        if key in R1_dict:
            del R1_dict[key]
        if key in R2_dict:
            del R2_dict[key]        
    print("number R1 reads after deleting empty and/or short read Pairs = " + str(len(R1_dict)))   
    print("number R2 reads after deleting empty and/or short read Pairs = " + str(len(R2_dict)))
    
    # where 1 pair is longer than cutoff and other is short or missing
    widow_R1 = []
    widow_R2 = []
    listDiffR1 = list(set(listComboR1) - set(listBoth))
    listDiffR2 = list(set(listComboR2) - set(listBoth))
    print("num R1 empty	or short where pair is long = "	+ str(len(listDiffR1)))
    print("num R2 empty	or short where pair is long = "	+ str(len(listDiffR2)))

    ## if R1 read is empty or short, then remove in R1 file, move R2 pair to widow, remove R2 pair
    for key in listDiffR1:
        if key in R2_dict:
            widow_R2.append(R2_dict[key])
            del R2_dict[key]
        if key in R1_dict:
            del R1_dict[key]
    ## if R2 read is empty or short, then remove in R2 file, move R1 pair to widow, remove R1 pair
    for key in listDiffR2:
        if key in R1_dict:
            widow_R1.append(R1_dict[key])
            del R1_dict[key]
        if key in R2_dict:
            del R2_dict[key]

    print("number widow R1 reads (where R2 is missing or short) = " + str(len(widow_R1)))   
    print("number widow R2 reads (where R1 is missing or short) = " + str(len(widow_R2)))  
    
    print("number R1 after deleting all empty/short reads = " + str(len(R1_dict)))   
    print("number R2 after deleting all empty/short reads = " + str(len(R2_dict)))   
       
    SeqIO.write(R1_dict.values(), args.widowDir + "/" + R1name + "_noEmpty.fastq", "fastq")
    SeqIO.write(R2_dict.values(), args.widowDir + "/" + R2name + "_noEmpty.fastq", "fastq")

    SeqIO.write(widow_R1, args.widowDir + "/" + R1name + "_widow.fastq", "fastq")
    SeqIO.write(widow_R2, args.widowDir + "/" + R2name + "_widow.fastq", "fastq")
    

if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()
