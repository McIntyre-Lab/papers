#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=gnsanujc
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/gnsanujc_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --array=1-500
##SBATCH --array=500
#SBATCH --qos=mcintyre-b

# run id_ujc.py by gene on all sample GTF

# dsan1: have 16075 unique genes in aln
# 16075/500 jobs -> 33 genes per job

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

# list of all genes (16075 rows)
DESIGN_FILE=$PROJ/design_files/list_dsan_data_gene.csv

## output dir
OUT=$PROJ/ujc_byGene_output
    mkdir -p $OUT

# cp gtf into slurm tempdir to be used for each job
if [ ! -f $SLURM_TMPDIR/dsan_data_2_dsan1_catTR.gtf ]; then
	cp $PROJ/gffcompare_by_TR/dsan_data_2_dsan1_catTR.gtf $SLURM_TMPDIR/dsan_data_2_dsan1_catTR.gtf
fi

INGTF=$SLURM_TMPDIR/dsan_data_2_dsan1_catTR.gtf

# set number of genes per job
GENES_PER_JOB=33
#GENES_PER_JOB=2

# Loop through above number of genes in one array job
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $GENES_PER_JOB + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $GENES_PER_JOB ))

for (( RUN=$START_NUM; RUN<=$END_NUM; RUN ++ )); do

	# RUN = GENE NUM (line number in design file)
	DESIGN=$(sed -n "${RUN}p" $DESIGN_FILE)
	IFS=',' read -ra ARRAY <<< "$DESIGN"

	GENE=${ARRAY[0]}

	if [ -z $GENE ]; then
		echo "Empty Gene."
	else
		date
		echo "ID UJC for $GENE in dsan_data"

		awk -v gene=${GENE} '$10~gene' \
			${INGTF} \
			> ${SLURM_TMPDIR}/dsan_data_${GENE}.gtf

		### ID UJC
		module purge
		export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

		python $UTIL/id_ujc_02ksb.py \
			-g $SLURM_TMPDIR/dsan_data_${GENE}.gtf \
			--keep-geneID \
			--track-source \
			-o $OUT

		date
		echo "$GENE! Wow! DoneDone!"

	fi
done



