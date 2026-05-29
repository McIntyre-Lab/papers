#!/bin/bash
#SBATCH --job-name=T1D_ujc
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/T1D_ujc.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=300gb
#SBATCH --time=96:00:00
#SBATCH --account=concannon
#SBATCH --qos=concannon-b

module purge

#Load modules
ml conda
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

##Create directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio

SCRIPTS=/blue/mcintyre/k.bankole/github/TranD/utilities

IND=${PROJ}/sam2gtf_all_runs/corrected_geneID_gtf_all_runs

OUTD=$PROJ/id_UJC
ROZ=${OUTD}/roz
    mkdir -p $OUTD
    mkdir -p ${ROZ}

##Cat all samples read gtfs together

cat ${IND}/*.gtf > ${ROZ}/allSamples.gtf

#
echo "Identifying UJCs from all Samples..."
## ID UJCs

python ${SCRIPTS}/id_ujc.py \
-g ${ROZ}/allSamples.gtf \
-o ${OUTD} \
-p allSamples

date

rm -r ${ROZ}
