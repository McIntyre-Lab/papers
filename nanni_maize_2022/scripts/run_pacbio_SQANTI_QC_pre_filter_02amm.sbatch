#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=SQ_pre
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/SQ_pre_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8

module purge
module load python/2.7.14
module load gcc/5.2.0
module load gmap/20180325
module load R/3.4.3
module load perl/5.24.1



### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
SCRIPTS=$PROJ/scripts/ConesaLab-sqanti

OUT=$PROJ/sqanti_qc_pre_filter
    mkdir -p $OUT

JUNCTIONS=$PROJ/star_dataPrep_junctions/SJout/
NAME="STAR"

## Abundance file list from tofu2 step --> cat contents of file (comma separated list) to variable
FL=$(cat $PROJ/tofu2_b73/abundance_file_list.csv)

## reference GTF file
REF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf
## reference fasta file
FA=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

echo  "$(date)

Junctions : ${JUNCTIONS}

abundance file :  $FL

Reference fasta : ${FA}

Reference GTF : ${REF}

"

### PRE-FILTER SQANTI

## Get chained PacBio isoforms GTF
##    (file ending is .gff but it is a gtf)
GTF=$PROJ/pacbio_analysis/merge_samples_b73/all_samples.chained.gff

## expression matrix
EXP=$PROJ/rsem_output/isoform_results/combined_expression_matrix.pre-filter.expected_count.txt

    echo "Running pre-filtered SQANTI...$(date)
    GTF used : ${GTF}"

    ## Make output directory if not already made
#    if [ ! -e $OUT/${NAME} ]; then
#        mkdir -p $OUT/${NAME}
#    fi
    ## Run pre-filtered SQANTI
    python ${SCRIPTS}/sqanti_qc.py \
        -g ${GTF} ${REF} ${FA} \
        -o SQANTI \
        -d $OUT \
        -c ${JUNCTIONS} \
        -fl ${FL} \
        -e ${EXP}
    echo "Done.

"
