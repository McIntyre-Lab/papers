#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=koppCatGTF
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/catGTF_%A_%a.out
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: cat noXLOC gtf into one big species GTF

## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
DIR=$PROJ/gffcompare_by_TR

SPECIES=dser
GENOME=dser1

echo $DIR
echo ${SPECIES}_*_TR*_2_${GENOME}_noXLOC.gtf

cat $DIR/${SPECIES}_*_TR*_2_${GENOME}_noXLOC.gtf > $DIR/${SPECIES}_data_2_${GENOME}_catTR.gtf

echo "Done Catting!"

echo ""

#echo "Sorting..."
# Sort by seqname, start, and transcriptID
#sort -t $'\t' -k1,1 -k4,4n -k10,10 $DIR/${SPECIES}_data_2_${GENOME}_catTR.gtf > $DIR/${SPECIES}_data_2_${GENOME}_catTR_sorted.gtf
