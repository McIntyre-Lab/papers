#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=mssngGn
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/mssngGn_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --account=mcintyre

module purge

## Description

export PATH=/ufgi-grimshawi/TB14/k.bankole/conda_env/TranD/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing/

for GENOME in dmel6 dsim2 dsan1 dyak2 dser1; do

#	echo "gene_id" > $PROJ/fiveSpecies_annotations/list_fiveSpecies_2_${GENOME}_genes_missing_from_event_analysis.csv
	> $PROJ/fiveSpecies_annotations/list_noHeader_fiveSpecies_2_${GENOME}_genes_missing_from_event_analysis.csv

	comm -23 <(tail -n +2 $PROJ/fiveSpecies_annotations/list_fiveSpecies_2_${GENOME}_ujc_uniq_gene_set.csv | sort) \
	<(tail -n +2 $PROJ/trand_1gtf_geneMode_fiveSpecies/fiveSpecies_2_${GENOME}_ujc_event_analysis_er.csv | awk -F, '{print $1}' | sort | uniq) \
	>> $PROJ/fiveSpecies_annotations/list_noHeader_fiveSpecies_2_${GENOME}_genes_missing_from_event_analysis.csv

	echo "$GENOME: "
 	echo $(tail -n +2  $PROJ/fiveSpecies_annotations/list_noHeader_fiveSpecies_2_${GENOME}_genes_missing_from_event_analysis.csv | wc -l)
done
