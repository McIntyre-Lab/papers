#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=parse_intpr
#SBATCH --output=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB/scripts/SLURM_LOGS/parse_interproscan_%A.out
#SBATCH --qos=mcintyre
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1

module load java/1.8.0_31 python perl

### Parse of InterProScan
### Combine all related parsed files into one
### 	(all files are still split)

### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
OUT=$PROJ/interproscan
SPECIES=maize

## Move all XML files into a directory by themselves
if [ ! -e ${OUT}/xml_output ]; then mkdir ${OUT}/xml_output; fi
mv ${OUT}/*.xml ${OUT}/xml_output

## Parse InterProScan Output
## Use lastparseInterproscanXml.py script provided by Fran
python /ufrc/mcintyre_conesa_transvar/isoAnnot/sorghum_maize_analysis/scripts/lastparseInterproscanXml.py \
    ${OUT}/xml_output ${OUT}/${SPECIES}_interproscan_out.tsv

## Cat all other analysis parsed files with interproscan parsed output
cat ${OUT}/${SPECIES}.??.*.parsed.out ${OUT}/${SPECIES}_interproscan_out.tsv \
	| sort -k1,1 >  ${OUT}/final.${SPECIES}_interproscan_out.tsv

