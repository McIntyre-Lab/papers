#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=SQANTI_QC_post
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/pacbio_sqanti_QC_post_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

module purge
module load python/2.7.14
module load gcc/5.2.0
module load gmap/20180325
module load R/3.4.3
module load perl/5.24.1


### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
SCRIPTS=$PROJ/scripts/ConesaLab-sqanti

OUT=$PROJ/sqanti_post_filter
    mkdir -p $OUT


JUNCTIONS=$PROJ/star_dataPrep_junctions/SJout/
NAME="STAR"

## Abundance file list from tofu2 step --> cat contents of file (comma separated list) to variable
FL=$(cat $PROJ/tofu2_b73/abundance_file_list.csv)

## ref GTF file
REF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf
## reference fasta file
FA=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa


    echo "$(date)

Junctions : ${JUNCTIONS}
Abundance file : ${FL}
Reference fasta is ${FA}
Reference GTF is ${REF}

"

### POST-FILTER SQANTI

    ## Make temporary directory
    ROZ=$OUT/roz_temp
        mkdir -p $ROZ

    ## Get filtered PacBio isoforms GTF -->> this is GTF from FILTERING step
    GTF=$PROJ/sqanti_filter/sqanti_filtered.gtf

    ## Combine all rsem expression files to make expression matrix -->> this is NEW expression matrix!!

    ## From Fran script : /ufrc/conesa/fpardopalacios/fly/rsem/mel/quantification/parse_expression.sh
    EXP=$PROJ/rsem_output_postFilter/isoform_results/combined_expression_matrix.post-filter.expected_count.txt

    echo "Running post-filtered SQANTI...$(date)
    GTF used : ${GTF}"


    ## Run post-filtered SQANTI
    python ${SCRIPTS}/sqanti_qc.py \
        -g ${GTF} ${REF} ${FA} \
        -o SQANTI \
        -d $OUT \
        -c ${JUNCTIONS} \
        -fl ${FL} \
        -e ${EXP}

    echo "Done."
    rm -r ${ROZ}

done
