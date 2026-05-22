#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=trand_keepIR
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/trand_keepIR_%A-%a.out
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --account=mcintyre
#SBATCH --array=0-4
##SBATCH --qos=mcintyre-b

module purge

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing/

IND=$PROJ/fiveSpecies_annotations

ARRAY=(dmel6 dsim2 dsan1 dyak2 dser1)

GENOME=${ARRAY[$SLURM_ARRAY_TASK_ID]}

OUTD=$PROJ/trand_1gtf_geneMode_fiveSpecies
	mkdir -p $OUTD

ORIG_GTF=$IND/fiveSpecies_2_${GENOME}_ujc.gtf
SUBS_GTF=$SLURM_TMPDIR/fiveSpecies_2_${GENOME}_ujc_keepIR.gtf

MSNG_GN_LST=$IND/list_noHeader_fiveSpecies_2_${GENOME}_genes_missing_from_event_analysis.csv

grep -f $MSNG_GN_LST $ORIG_GTF > $SUBS_GTF

date

NAME=fiveSpecies_2_${GENOME}_ujc_keepIR

trand \
$SUBS_GTF \
-o $OUTD \
-e gene \
-x ${NAME} \
-f \
-k \
-n $SLURM_CPUS_ON_NODE \
--verbose \
-l $OUTD/${NAME}.log

date

