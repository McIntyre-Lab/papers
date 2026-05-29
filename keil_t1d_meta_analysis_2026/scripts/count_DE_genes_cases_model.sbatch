#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=DE_case.%j
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/DE_case.%A_%a
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-5


###Set directories and files
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IND=${PROJ}/model_DE_frag_case_CD4_CD8/
OUTD=${IND}/DE_gene_counts
    mkdir -p ${OUTD}
ROZ=$OUTD/roz
    mkdir -p ${ROZ}


## design has header!!
DESIGN_FILE=$PROJ/design_files/design_effect_controls.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
EFFECT=${ARRAY[1]}

T1D=${PROJ}/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt
IMMUNE=${PROJ}/candidate_gene_lists/immunobase_hg19_ensembl_entrez.csv

SIG=0.001
PROP=0.2

module load python/3.6

python $PROJ/scripts/count_DE_genes_from_frag_model.py \
        --t1d ${T1D} \
        --immune ${IMMUNE} \
        --in_file ${IND}/de_t1_frag_case.csv \
        --effect ${FLAG} \
        --pvalue ${SIG} \
        --prop ${PROP} \
        --out_graph ${OUTD}/prop_gene_DE_${EFFECT}_effect_sig_${SIG}_cut_${PROP}.png \
        --out_csv ${OUTD}/gene_DE_cnts_${EFFECT}_effect_sig_${SIG}_cut_${PROP}.csv
