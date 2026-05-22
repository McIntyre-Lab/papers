#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=5s_eres
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/5s_eres_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --account=mcintyre
#SBATCH --array=0-4

module purge

export PATH=/ufgi-grimshawi/TB14/k.bankole/conda_env/TranD/bin:$PATH
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

GENOME_ARR=(dmel6 dsim2 dsan1 dyak2 dser1)

GENOME=${GENOME_ARR[$SLURM_ARRAY_TASK_ID]}

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing

IND=$PROJ/trand_1gtf_geneMode_fiveSpecies
OUTD=$PROJ/fiveSpecies_annotations

cp $IND/fiveSpecies_2_${GENOME}_ujc_event_analysis_er.csv $SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_event_analysis_er.csv
tail -n +2 $IND/fiveSpecies_2_${GENOME}_ujc_keepIR_event_analysis_er.csv >> $SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_event_analysis_er.csv

cp $IND/fiveSpecies_2_${GENOME}_ujc_event_analysis_ef.csv $SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_event_analysis_ef.csv
tail -n +2 $IND/fiveSpecies_2_${GENOME}_ujc_keepIR_event_analysis_ef.csv >> $SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_event_analysis_ef.csv

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

python $UTIL/convert_event_analysis_2_gene_gtf_01ksb.py \
	-e $SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_event_analysis_er.csv \
	-o $OUTD/fiveSpecies_2_${GENOME}_ujc_er.gtf

python $UTIL/convert_exon_segment_2_es_gtf_01ksb.py \
	-e $SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_event_analysis_ef.csv \
	-o $OUTD/fiveSpecies_2_${GENOME}_ujc_es.gtf

