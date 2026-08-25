#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=gnmelujc
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/gnmelujc_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --array=1-500
##SBATCH --qos=mcintyre-b

# run id_ujc.py by gene on all sample GTF
# if rerun, DELETE OLD OUTPUT

# dmel6: have 16712 unique genes in aln
# 16712/500 jobs -> 34 genes per job

## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

# list of all genes (16712 rows)
DESIGN_FILE=$PROJ/design_files/list_mel_data_gene.csv

## output dir
OUT=$PROJ/ujc_byGene_output
    mkdir -p $OUT

# cp gtf into slurm tempdir to be used for each job

if [ ! -f $SLURM_TMPDIR/mel_data_2_dmel6_catTR.gtf ]; then
	cp $PROJ/gffcompare_by_TR/mel_data_2_dmel6_catTR.gtf $SLURM_TMPDIR/mel_data_2_dmel6_catTR.gtf
fi

INGTF=$SLURM_TMPDIR/mel_data_2_dmel6_catTR.gtf

# set number of genes per job
GENES_PER_JOB=34
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
		echo "ID UJC for $GENE in dmel_data"

		awk -v gene=${GENE} '$10~gene' \
			${INGTF} \
			> ${SLURM_TMPDIR}/dmel_data_${GENE}.gtf

		### ID UJC
		module purge
		export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

		python $UTIL/id_ujc_02ksb.py \
			-g $SLURM_TMPDIR/dmel_data_${GENE}.gtf \
			--keep-geneID \
			--track-source \
			-o $OUT

		date
		echo "$GENE! Wow! DoneDone!"
		rm $SLURM_TMPDIR/dmel_data_${GENE}.gtf
	fi
done



