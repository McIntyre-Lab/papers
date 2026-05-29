import pandas as pd
import argparse
import os
import sys
import csv
from Bio import SeqIO

def getoptions():
    parser = argparse.ArgumentParser(description='Extract counts information from lima outputs')
    parser.add_argument("-l", "--list", dest="list", required=True, help="file that contains a list of transcript names")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output file name with full path")
    parser.add_argument("-i", "--input", dest="input", required=True, help="input file name with full path")
    parser.add_argument("-t", "--type", dest="type", required=True, default= "sam", help="extract information based on a gene/transcript list from assigned file type. The type can be sam, fasta and classification")
    args = parser.parse_args()
    return(args)


def parseSamHeader(sam):
    counts = 0
    with open(sam, "r") as f:
        for line in f:
            if line.startswith("@"):
                counts+=1
            else:
                return counts


def extractFromSam(list, sam, out):
    n_header = parseSamHeader(sam)
    transcript_list = pd.read_table(list, header = None, names = ["name"])
    samfile = pd.read_csv(sam, sep = "\n", header = None, names = ["name"])
    header_df = samfile.iloc[range(0, n_header),:]
    namecol = samfile.name.str.split("\t").str[0]
    newSAM = samfile[namecol.isin(transcript_list.name)]
    newSAM = header_df.append(newSAM)
    newSAM.to_csv(out, sep="\n", index = False, header = False, quoting=csv.QUOTE_NONE, quotechar = "", escapechar = '\n')
    return
                    


def extractFromFasta(list, fa, out):
    transcript_list = pd.read_table(list, header = None, names = ["name"])
    identifiers = [record.id for record in SeqIO.parse(fa, "fasta")]
    sequences = [str(record.seq) for record in SeqIO.parse(fa, "fasta")]
    df = pd.DataFrame({"name": identifiers, "sequence": sequences})
    new_df = df[df.name.isin(transcript_list.name)]
    new_df["name"] = ">" + new_df["name"]
    new_df.to_csv(out, sep="\n", index = False, header = False)
    return

def extractFromClass(list, cls, out):
    transcript_list = pd.read_table(list, header = None, names = ["name"])
    cls_df = pd.read_table(cls, header = 0, sep = "\t", dtype=str)
    new_cls = cls_df[cls_df.isoform.isin(transcript_list.name)]
    new_cls.columns = cls_df.columns
    new_cls.to_csv(out, sep = "\t", na_rep = "NA", index = False)
    return


def main():
    args = getoptions()
    if args.type == 'sam':
        extractFromSam(args.list, args.input, args.out)

    elif args.type == 'fasta':
        extractFromFasta(args.list, args.input, args.out)

    elif args.type == 'class':
        extractFromClass(args.list, args.input, args.out)
    else:
        sys.stderr.write("Invalid type assigned! Only sam, fasta and class are valid")
    return

if __name__ == "__main__":
    main()
