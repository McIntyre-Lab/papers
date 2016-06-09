import os
import pandas as pd
import mclib.vcf2 as mcvcf

mclab = os.getenv('MCLAB')
sandbox = '/home/jfear/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files'

## File Names ##
design = os.path.join(mclab, 'cegs_ase_paper/design_files/CEGS_list_68_lines.txt')
output = os.path.join(mclab, 'cegs_ase_paper/pipeline_output/polymorphisms_counts.csv')

## Prepare result dataset ##
# Create a list of lines
with open(design, 'r') as DS:
    lines = [x.strip() for x in DS]

# Create a results dataset
results = pd.DataFrame(index=lines, columns=['masked', 'snps', 'indels', 'poly', 'total'], dtype=int)
results.index.name = 'line'

## Count polymorphisms ##


def countMasked(line, results):
    """ Function to count the number of row in the perm mask bed file """
    global sandbox
    fname = os.path.join(sandbox, 'pos_to_permMask_w11182{0}.bed'.format(line))
    with open(fname, 'r') as BED:
        cnt = 0
        for row in BED:
            cnt += 1
        results.ix[line, 'masked'] = cnt


def countVcf(line, results):
    """ Count the number of snps and indels in the vcf file """
    global sandbox
    fname = os.path.join(sandbox, '{0}_w11182{0}_UPD.vcf.gz'.format(line))
    vcf = mcvcf.Vcf(fname)
    scnt = 0
    icnt = 0
    for row in vcf.vcf_reader:
        if row.is_snp:
            scnt += 1
        elif row.is_indel:
            icnt += 1
    results.ix[line, 'snps'] = scnt
    results.ix[line, 'indels'] = icnt


## Iterate over lines and fill in results table with counts ##
# Here I count the number of locations that are permanently masked due to
# overlapping polymorphisms. I also count the number of SNPs and INDELs.

for line in lines:
    countMasked(line, results)
    countVcf(line, results)

results['poly'] = results['snps'] + results['indels']
results['total'] = results['poly'] + results['masked']

## Write out results ##
results.to_csv(output)
