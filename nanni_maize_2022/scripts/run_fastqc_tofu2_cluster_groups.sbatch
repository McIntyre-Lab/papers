#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=fastqc
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/fastqc_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --array=2-11

module load fastqc/0.11.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth

#for i in b73 mo17_yan mo17_cau
for REF in b73
do
    ## Set directories
    TOFU=$PROJ/compare_b73_2_mo17/tofu2_${REF}
    IGNORED=$PROJ/compare_b73_2_mo17/uncollapsed_${REF}_IsoSeq_reads
    OUTD=$PROJ/compare_b73_2_mo17/fastqc_${REF}_tofu2_groups
        mkdir -p ${OUTD}

    ## Design file
    DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    GENO=${ARRAY[1]}
    TRT=${ARRAY[2]}    
    ID=${ARRAY[3]}

    SAMPLE=${ID}_${GENO}_${TRT}

    ## Loop over the 3 groups of clusters following tofu2 collapse
    for TYPE in low_cov low_ident pass_tofu_filter; do
        ## Get FASTQ file
        if [[ ${TYPE} == "pass_tofu_filter" ]]; then
            FQ=${TOFU}/${SAMPLE}/${SAMPLE}.collapsed.rep.fq
        else
            FQ=${IGNORED}/${SAMPLE}_unCollapsed_${TYPE}_${REF}_ref.fq
        fi

        ## Run FastQC
        fastqc \
            --outdir ${OUTD} \
            ${FQ}

    done
done
