#!/usr/bin/evn python
""" A set of functions for calculating haplotype frequencies.

Table 5 in the manuscript has a summary of haplotype frequencies. This set of
functions are used by ipython notebooks to create the counts for this table.

"""

import vcf as pyvcf
import bed as mcbed
import pandas as pd


def convertLine(x):
    """ Converts CEGS genotype names to match vcf.

    This function converts genotype names from this format:
        r###
        w###

    To the format used by the vcf files:
        Raleigh_###
        Winters_###

    :param str x: The genotype name to be converted.

    :rtype: str
    :returns: The converted string.

    """
    if x.startswith('r'):
        x = x.replace('r', 'Raleigh_')
    elif x.startswith('w'):
        x = x.replace('w', 'Winters_')
    return x.strip()


def cegsAlleleCnt(bedName, vcfName, genotypes):
    """ Creates a table of haplotypes for genes of interest.

    Pulls out a list of genes and gentoypes from a large
    vcf file and counts the number of haplotypes for the
    genes of interest.

    :param str bedName: BED formatted file with genes of interest

    :param str vcfName: VCF file to parse, needs to be compressed using
        bgzip and indexed using tabix. Both of these tools are part of samtools.

    :param list genotypes: List of genotypes of interest

    :rtype: list of tuple
    :returns: A list of tuples with (gene, haplotype count)

    """
    # Read in a bed file
    bed = mcbed.Bed(bedName)

    # Connect to vcf file. Needs to be compressed with bgzip and indexed with tabix
    vcfReader = pyvcf.VCFReader(filename=vcfName, compressed=True)

    outList = list()
    # Iterate over bed file
    for row in bed:
        gene = row[3]
        # Fetch SNPs for current region
        for snp in vcfReader.fetch(row[0], row[1], row[2]):
            refBase = snp.REF
            pos = snp.POS

            # Get base for each genotype of interest
            for geno in genotypes:
                # There are not snp calls for Raleigh_286 so skip
                if geno != 'Raleigh_286':
                    curGeno = snp.genotype(geno)

                    if curGeno.called:
                        curBase = curGeno.gt_bases[0]
                    else:
                        # If SNP call missing then set to ref base
                        curBase = refBase
                    outList.append((gene, pos, refBase, geno, curBase))

    # Make a data frame
    df = pd.DataFrame(outList, columns=['gene', 'pos', 'ref', 'genotype', 'base'])
    df.set_index(['gene', 'pos', 'ref', 'genotype'], inplace=True)

    # Make a side-by-side
    # genotype as column headers and base as value
    dfWide = df.unstack(level='genotype')

    # For each gene figure out unique haplotypes
    grp = dfWide.groupby(level='gene')

    out = list()
    for name, g in grp:
        uniq = set([tuple(x) for x in g.T.values])
        out.append([name, len(uniq)])

    return out


