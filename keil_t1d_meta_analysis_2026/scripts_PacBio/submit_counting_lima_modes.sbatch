#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=count_demux_modes
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/count_demux_modes_%A-%a.out
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-4

module load samtools/1.12

### Set Project Directory
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
OUTD=${PROJ}/check_lima_modes
    mkdir -p ${OUTD}

## design has header!!
DESIGN_FILE=$PROJ/design_files/design_lima_modes.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

MODE=${ARRAY[0]}
PROC_DIR=${ARRAY[1]}
REFINE_DIR=${ARRAY[2]}



for RUN in first_run second_run third_run fourth_run fifth_run
do

if [[ ${MODE} == "symmetric" ]] ; then
    echo "mode is ${MODE}"


    num_ccs=$(samtools view -c /blue/concannon/share/T1D_R01_controls_PacBio/isoseq3_library_processing/${RUN}/*.ccs.bam)

    echo "starting number of ccs reads is ${num_ccs}" > $OUTD/starting_ccs_${RUN}.txt


    for DEMUX_FILE in ${PROJ}/${PROC_DIR}/${RUN}/demultiplexed/*.bam 
    do 
        lima_reads=$(samtools view -c ${DEMUX_FILE})
       	demux_filename=$(basename ${DEMUX_FILE} .bam)

       	echo "${demux_filename}:  ${lima_reads}" >> ${OUTD}/lima_counts_${MODE}_${RUN}.txt
    done
    
    for REF_FILE in ${PROJ}/${REFINE_DIR}/${RUN}/*.bam
    do
      	refine_reads=$(samtools view -c ${REF_FILE})
        ref_filename=$(basename ${REF_FILE} .bam)

        echo "${ref_filename}:  ${refine_reads}" >> ${OUTD}/isoseq3_refine_counts_${MODE}_${RUN}.txt
    done

   

else
    echo "mode is ${MODE}"
    for DEMUX_FILE in ${PROJ}/${PROC_DIR}/${RUN}/*.bam
    do
        lima_reads=$(samtools view -c ${DEMUX_FILE})
        demux_filename=$(basename ${DEMUX_FILE} .bam)

        echo "${demux_filename}:  ${lima_reads}" >> ${OUTD}/lima_counts_${MODE}_${RUN}.txt 
    done

    for REF_FILE in ${PROJ}/${REFINE_DIR}/${RUN}/*.bam
    do
      	refine_reads=$(samtools view -c ${REF_FILE})
        ref_filename=$(basename ${REF_FILE} .bam)

        echo "${ref_filename}:  ${refine_reads}" >> ${OUTD}/isoseq3_refine_counts_${MODE}_${RUN}.txt
    done

fi

done
