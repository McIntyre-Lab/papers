#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=compare_miso_bwa
#SBATCH --partition=hpg1-compute
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.compare_miso_bwa.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --array=2-6
#

#6
## Modules

module load miso/0.5.3
module load samtools
module load bedtools
module load python/2.7.6

## Set directories
    PROJ=/ufrc/mcintyre/share/conesa_isoform_check
    ORIG=$PROJ/original_data

    OUTPUT=$PROJ/bwa_mem_aln_pe
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

     ## Design file
    DESIGN_FILE=$PROJ/design_files/miso_compare_design.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE1=${ARRAY[0]}
    SAMPLE2=${ARRAY[1]}

### Pairwise comparison of MISO output

    OUTPUT=$PROJ/comparisons
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

     echo "Running MISO"

     SAMP1DATA=$PROJ/miso_output_SE_${SAMPLE1}
     SAMP2DATA=$PROJ/miso_output_SE_${SAMPLE2}


compare_miso --compare-samples ${SAMP1DATA} ${SAMP2DATA} $OUTPUT

