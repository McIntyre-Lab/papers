#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=full_pass
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/full_pass_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=1
#SBATCH --account=concannon
#SBATCH --qos=concannon-b

module purge
module load R/4.1 samtools/1.10

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
ROZ=${PROJ}/roz
IND=${PROJ}/isoseq3_library_processing

module load R/4.1 samtools/1.10

## Get distribution of full passes in CCS reads

## Each ccs.bam file has a “np:i:” flag in the 14th column that lists
##     the number of full passes used for the read

## the max is 60, because even if there are more than 60, ccs will 
##     just use the longest 60 subreads to generate the ccs read
##     (See https://ccs.how/faq/accuracy-vs-passes.html#is-there-an-upper-limit-on-number-of-passes-used)

## Set directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IND=$PROJ/isoseq3_library_processing
OUT=${IND}/num_full_pass_distribution.txt

#rm ${OUT}

for RUN in first second third fourth fifth;
do
   if [[ ${RUN} == 'first' ]]
   then

   BAM=${IND}/*ccs.bam

    ## Get distribution for first run
    echo "T1D_R01 ${RUN} run:" >> ${OUT}
    samtools view ${BAM} | cut -f 14 | sed s/"np:i:"/""/g \
            | Rscript -e 'summary (as.numeric (readLines ("stdin")))' \
            >> ${OUT}
        echo "
" >> ${OUT}

   else
   BAM=${IND}/${RUN}_run/*ccs.bam

    ## Get distribution for other runs
    echo "T1D_R01 ${RUN} run:" >> ${OUT}
    samtools view ${BAM} | cut -f 14 | sed s/"np:i:"/""/g \
            | Rscript -e 'summary (as.numeric (readLines ("stdin")))' \
            >> ${OUT}
        echo "
" >> ${OUT}
    fi
done
