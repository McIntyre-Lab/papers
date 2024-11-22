#!/bin/bash

## Description: Remove gene features from NCBI genomic annotation and subset dmel annotation for transpliced gene

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
SCRIPTS=$PROJ/scripts
REF=/blue/mcintyre/share/references

# Remove gene features from NCBI annotations
awk '$3 != "gene" ' $REF/dser1.1/GCF_002093755.2/genomic.gtf > $REF/dser1.1/GCF_002093755.2/genomic_gtf_no_gene_features.gtf
awk '$3 != "gene" ' $REF/dsan_Prin_1.1/GCF_016746245.2/genomic.gtf > $REF/dsan_Prin_1.1/GCF_016746245.2/genomic_gtf_no_gene_features.gtf
awk '$3 != "gene" ' $REF/dyak_Prin_Tai18E2_2.1/GCF_016746365.2/genomic.gtf > $REF/dyak_Prin_Tai18E2_2.1/GCF_016746365.2/genomic_gtf_no_gene_features.gtf

# Subset dmel650 for one transpliced gene
echo "FBgn0002781" > $PROJ/gene_exclusion.txt

python /blue/mcintyre/k.bankole/github/TranD/utilities/subset_gtf.py \
	-g /blue/mcintyre/share/references/dmel_fb650/dmel-all-r6.50.gtf \
	-t gene_id \
	-e $PROJ/gene_exclusion.txt \
	-o /blue/mcintyre/share/references/dmel_fb650/dmel-all-r6.50_subset_trans_spliced_gene.gtf

rm $PROJ/gene_exclusion.txt
