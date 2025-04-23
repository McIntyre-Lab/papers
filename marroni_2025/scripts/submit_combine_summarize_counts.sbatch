#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=prep
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/prep.%A-%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb

## combine ase count tables
## summarize data

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
SCRIPTS=$PROJ/scripts
SAMC=$PROJ/ase_counts_2_upd_w1118_and_2_upd_ril_02amm

SSUM=$PROJ/ase_counts_summed_techReps
    mkdir -p $SSUM
FILT=$PROJ/ase_counts_summarized
    mkdir -p $FILT

BEDFILE=/blue/mcintyre/share/references/dmel_fb643/dmel643_gene_region_plus_minus_100.bed
DESIGN=$PROJ/design_files/df_all_runs_samcompare_nfn03.csv

    # Load modules
    module load python/3.6 htslib

    mkdir -p ${SSUM}

    SIM=False
    APN=1

    python3 ${SCRIPTS}/combine_cnt_tables_13amm.py \
        -design ${DESIGN} \
        -sim ${SIM} \
        --bed ${BEDFILE} \
        --path ${SAMC} \
        --designdir $PROJ \
        --out ${SSUM}


    ## design2 is created in above script
    DESIGN2=$PROJ/df_ase_samcomp_summed.csv

    mkdir -p ${FILT}

    python3 ${SCRIPTS}/summarize_sam_compare_cnts_table_1cond_and_output_APN_06amm.py \
        --output ${FILT} \
        --design ${DESIGN2} \
        --parent1 G1 \
        --parent2 G2 \
        --sampleCol sample \
        --sampleIDCol sampleID \
        --sam-compare-dir ${SSUM} \
        --apn ${APN} \
