#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=catGTFK
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/catGTF_%A-%a.out
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

# Cat all samples per species -> one GTF.

module purge

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

OUT=$PROJ/gffcompare_all_sample_per_species
	mkdir -p $OUT

SPECIES=dser
GENOME=dser1

# Cat all samples to one GTF.

date
echo "Catting all samples..."

cat $OUT/${SPECIES}_*_2_${GENOME}.gtf > $OUT/${SPECIES}_data_2_${GENOME}.gtf

date

