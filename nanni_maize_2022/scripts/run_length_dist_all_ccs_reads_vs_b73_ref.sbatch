#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=length_ccs
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/length_ccs_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Get lengths of all ccs reads (after trimming)

## Plot distribution of ccs read and reference transcript lengths

module load gffread/0.12.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/check_isoseq3_lost_genes/bam2fq
OUTD=$PROJ/compare_b73_2_mo17/ccs_read_length
    mkdir -p ${OUTD}
ROZ=${OUTD}/roz_ccs_length
    mkdir -p ${ROZ}

## Get reference transcript fasta
REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
REF_GTF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf

gffread -w ${ROZ}/ref_b73.fa \
    -g ${REF} ${REF_GTF}

## Set up env to have seaborn v0.9.0 for python code to run properly
export PATH=/blue/mcintyre/adalena.nanni/conda_envs/plot_seaborn_0.9.0/bin:$PATH

## Concatenate samples within each genotype
for GENO in b73 c123 hp301 mo17 nc338; do
    echo "Evaluate ${GENO} ccs reads"
    cat ${IND}/*_${GENO}*.ccs.trimmed.primer_5p--primer_3p.fq \
        > ${ROZ}/${GENO}_ccs.fq

    ## Too many reads in FASTQ to efficiently process in python script
    ## Getting read lengths using awk
    awk 'BEGIN{print "length";} NR%4==2{print length($1)}' \
        ${ROZ}/${GENO}_ccs.fq \
        > ${ROZ}/${GENO}_ccs.csv
    rm ${ROZ}/${GENO}_ccs.fq

    ## Plot distribution of ccs sequence lengths vs. reference lengths
    echo "...plot"
    if [[ ! -e ${OUTD}/ref_b73_seq_length.csv ]]; then
        python $PROJ/scripts/plot_seq_length_distribution.py \
            -l ${ROZ}/${GENO}_ccs.csv \
            -fa ${ROZ}/ref_b73.fa \
            -o ${OUTD} \
            -p ${GENO}_vs_b73_ref
    else
        python $PROJ/scripts/plot_seq_length_distribution.py \
            -l ${ROZ}/${GENO}_ccs.csv \
            -l ${ROZ}/ref_b73_seq_length.csv \
            -o ${OUTD} \
            -p ${GENO}_vs_b73_ref
    fi
done

rm -r ${ROZ}
