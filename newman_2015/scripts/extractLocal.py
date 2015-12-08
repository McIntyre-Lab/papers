#!/usr/bin/env python

#	DESCRIPTION: This program gets co-ordinates for pulling SNPs local to a gene for cis-eQTL.
#	The output is a TSV file with the list of 5' start, transcriptional start site (TSS)
#       region, stop position region for transcripts, and 3' stop. The 5' start and 3' stop positions
#       define the region considered "local" for a gene
#
#	
#   AUTHOR: Jeremy Newman

# Build-in packages
import argparse  # Command line use

# Add-on packages
import gffutils
import itertools
import numpy as np

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--input", dest="gffFile", action='store', required=True, help="Input database")
    parser.add_argument("--output", dest="outputFile", action='store', required=True, help="Output file name")
    parser.add_argument("--5prime", dest="upstreamSize", action='store', type=int, required=True, help="Size of 5-prime region in bp")
    parser.add_argument("--3prime", dest="downstreamSize", action='store', type=int, required=True, help="Size of 3-prime region in bp")
    parser.add_argument("--gene_list", dest="geneList", action='store', required=True, help="Genes to extract")

    args = parser.parse_args()
    return args

def main():
    # Get database file
    db_fn = args.gffFile
    db = gffutils.FeatureDB(db_fn)
    genes = db.features_of_type('gene', order_by='start')
    
    genelist=np.genfromtxt(args.geneList, delimiter=',', dtype=None)
    genelist_tpose=np.transpose(genelist)
    
    # Open output file
    with open(args.outputFile, 'wb') as outputFile:
        # Write column headings
        
        for gene in genes:
            for i in range(0,len(genelist_tpose)):
                if gene.id == genelist_tpose[i]:
                    gene_start = gene.start
                    gene_stop = gene.end 
                    genename = gene.id
                    genechrom = gene.chrom

                    # Upstream and downstram positions indicated by command line arguements
                    totalstart=gene_start - args.upstreamSize
                    totalstop=gene_stop + args.downstreamSize
                    bedArray = [genechrom, totalstart, totalstop, genename]
                    outputFile.write("\t".join(str(i) for i in bedArray) + "\n")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

