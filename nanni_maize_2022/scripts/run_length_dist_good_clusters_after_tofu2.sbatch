#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=length_good
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/length_good_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Get lengths of good clusters (pass tofu2 filtering)

## Plot distribution of cluster lengths and count how many clusters are
##     <500nt, >= 500nt and <1000nt, and >=1000nt

module load python/3.8
module load seqtk

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/compare_b73_2_mo17/blast_low_ident_cov_eval
    mkdir -p ${OUTD}
ROZ=${OUTD}/roz_good_cluster
    mkdir -p ${ROZ}

#for REF in b73 mo17_yan mo17_cau; do
for REF in b73; do

    ## Set input directory
    IND=$PROJ/compare_b73_2_mo17/tofu2_${REF}

    for FILE in $(ls ${IND}/*/*.collapsed.rep.fq); do
        SAMPLE=$(basename ${FILE} .collapsed.rep.fq)
        echo "Evaluate ${SAMPLE} good clusters mapped to ${REF}"

        seqtk seq -A ${FILE} > ${ROZ}/${SAMPLE}_2_${REF}.fa

        ## Plot distribution of sequence lengths
        python $PROJ/scripts/fasta_seq_length_distribution.py \
            -fa ${ROZ}/${SAMPLE}_2_${REF}.fa \
            -o ${OUTD}/${SAMPLE}_pass_tofu_filter_${REF}_ref

    done

done

rm -r ${ROZ}
