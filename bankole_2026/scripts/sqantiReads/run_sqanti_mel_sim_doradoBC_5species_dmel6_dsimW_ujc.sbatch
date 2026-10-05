#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sqR
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/sqR_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
##SBATCH --array=4
#SBATCH --array=5-109


## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
PYCHOP=${PROJ}/dorado_basecalling/pychopper
ALN=${PROJ}/aln_mel_sim_dorReads_2_dmel6_dsimW

SQANTI=${PROJ}/sqanti3QC_mel_sim_dorReads_2_5species_dmel6_dsimW_ujc
    mkdir -p $SQANTI

CNTS=${PROJ}/gene_fsm_ism_counts_from_sqanti_classFile
    mkdir -p $CNTS

## Get design file information (2 - 109)
DESIGN_FILE=$PROJ/design_files/sample_bc_design_w_origDataPath_04amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
BC=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[3]}
TECH=${ARRAY[4]}
DATE=${ARRAY[5]}
ONT_DIR=${ARRAY[6]}

SAMPLEID=${SAMPLE}_TR${TR}

ROZ=$PROJ/roz_sq_${SAMPLEID}
    mkdir -p $ROZ


if [[ $SAMPLEID == mel* ]]; then
    SAM=${ALN}/${SAMPLEID}_fl_rsc_2_dmel6.sam
    BED=${ALN}/${SAMPLEID}_fl_rsc_2_dmel6.bed
    GTF=${ALN}/${SAMPLEID}_fl_rsc_2_dmel6.gtf
    BAM=${ALN}/${SAMPLEID}_fl_rsc_2_dmel6.bam
    BAMS=${ALN}/${SAMPLEID}_fl_rsc_2_dmel6_sort.bam
    TAG=dmel6
else
    SAM=${ALN}/${SAMPLEID}_fl_rsc_2_dsimW.sam
    BED=${ALN}/${SAMPLEID}_fl_rsc_2_dsimW.bed
    GTF=${ALN}/${SAMPLEID}_fl_rsc_2_dsimW.gtf
    BAM=${ALN}/${SAMPLEID}_fl_rsc_2_dsimW.bam
    BAMS=${ALN}/${SAMPLEID}_fl_rsc_2_dsimW_sort.bam
    TAG=dsimW
fi

    module purge
    module load samtools/1.10 bedtools/2.29.2

    ## SAM to BAM

    samtools view -b -F 2048 ${SAM} > ${BAM}

    ## BAM to BED
    bedtools bamtobed -split -i ${BAM} > ${BED}

    ## BED to GTF
    SOURCE=${SAMPLEID}_fl_rsc
    awk -v source=${SOURCE} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' $BED > $GTF

    ## Remove intermediate files
    rm ${BED} ${BAM}


    ## set up refs to use for sqanti3qc
    if [[ $TAG == dmel6 ]]; then
        REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf
        REF_FA=/blue/mcintyre/share/references/dmel_fb650/dmel-all-chromosome-r6.50.fasta
    elif [[ $TAG == dsimW ]]; then
        REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dsimW_ujc.gtf
        REF_FA=/blue/mcintyre/share/references/dsim_WXD1_ASM438218v1/GCA_004382185.1_ASM438218v1_genomic.fasta
    else
        echo "ref setup for $SAMPLEID is not right!"
    fi

    module purge
    module load sqanti3/4.2

    sqanti3_qc.py \
        --skipORF \
        --report pdf \
        --force_id_ignore \
        --dir ${ROZ} \
        -o ${SAMPLEID}_fl_rsc_${TAG} \
        ${GTF} \
        ${REF_GTF} \
        ${REF_FA}

    ## Move all files to output directory
    mv ${ROZ}/*.* ${SQANTI}/.

    rm -r $ROZ

    ## transcript counts from sqanti class file
    module load python/3.10

    python $PROJ/scripts/count_FSM_ISM_from_classification.py \
        --classification ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_classification.txt \
        --out $CNTS/${SAMPLEID}_fl_rsc_${TAG}_FSM_ISM_Cnts.txt


    # add sampleID column to output
    awk -v OFS=',' -v sampleid=$SAMPLEID '{if (NR==1) print "sampleID", $0; else print sampleid, $0}' $CNTS/${SAMPLEID}_fl_rsc_${TAG}_FSM_ISM_Cnts.txt \
        > $CNTS/${SAMPLEID}_fl_rsc_${TAG}_FSM_ISM_Cnts.csv

    rm $CNTS/${SAMPLEID}_fl_rsc_${TAG}_FSM_ISM_Cnts.txt


    ## gene counts
    python $PROJ/scripts/gene_counts_from_classsification.py \
        --classification ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_classification.txt \
        --out $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt


    # add sampleID column to output
    awk -v OFS=',' -v sampleid=$SAMPLEID '{if (NR==1) print "sampleID", $0; else print sampleid, $0}' $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt \
        > $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.csv

    rm $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt

