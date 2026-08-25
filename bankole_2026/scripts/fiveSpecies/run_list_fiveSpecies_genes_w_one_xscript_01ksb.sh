#!/bin/bash

## List genes with one transcript

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

for GENOME in dmel6 dsim2 dsan1 dyak2 dser1; do
	python $UTIL/find_genes_with_one_transcript_01ksb.py \
		-g $PROJ/fiveSpecies_annotations/fiveSpecies_2_${GENOME}_ujc.gtf \
		-o $PROJ/fiveSpecies_annotations
done
