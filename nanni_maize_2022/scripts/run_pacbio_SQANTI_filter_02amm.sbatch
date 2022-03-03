#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=SQANTI_filter
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/pacbio_sqanti_filter_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=1
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
SCRIPTS=/ufrc/mcintyre/share/etoh_srna/scripts/ConesaLab-sqanti
##SCRIPTS=$PROJ/scripts/ConesaLab-sqanti

OUT=$PROJ/sqanti_filter
    mkdir -p $OUT

JUNCTIONS=/ufrc/mcintyre/share/maize_ainsworth/star_dataPrep_junctions/SJout/

echo "date

Junctions : ${JUNCTIONS}

"

## Make output directory if not already made 

    if [ ! -e $OUT/filter_out ]; then
        mkdir -p $OUT/filter_out
    fi
     
    ## Get PacBio GTF --> this is GTF from PRE-FILTERING step
    GTF=$PROJ/sqanti_qc_pre_filter/all_samples.chained_corrected.gtf

    ## Use SQANTI QC output in SQANTI filter to filter out probably artifacts
    echo "Filter transcriptome... $(date)
    filtering with GTF : $GTF"    
    python ${SCRIPTS}/sqanti_filter.py \
        /ufrc/mcintyre/share/maize_ainsworth/sqanti_qc_pre_filter/SQANTI_classification.txt \
        -d $OUT/filter_out
    echo "Done.
"
    CUR=$OUT/filter_out/*_curatedTranscriptome.txt
    awk -F '"' 'NR==FNR{ids[$1]=1;next} ids[$2]' ${CUR} ${GTF} \
        > $OUT/sqanti_filtered.gtf

