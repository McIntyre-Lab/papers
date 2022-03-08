#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
SCRIPTS=~/mclab/SHARE/McIntyre_Lab/useful_mclab_info/scripts/python
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
DESIGN=$PROJ/design_files
INM=$PROJ/cvr_cnts_gene_mo17_cau
INB=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file
OUTD=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/genetics_2022/B73_Mo17_CAU_cvr_cnt_mappedRead_scatter_12604
    mkdir -p ${OUTD}

## For the 12604 assembled transcriptome
## Identify the B73-Mo17 syntenic genes
## Plot scatter plot comparing the mapped reads when mapped to B73 vs. Mo17 CAU
## Make 2 plots per sample (genotype-treatment):
##     1) with techreps plotted separately (1 point per techrep-gene pair),
##     2) with techreps combined (1 point per gene)

## Get list of 12604 genes
GENE=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/resubmission_2021/supplement/Supplementary_File_5.csv

## Get B73-Mo17CAU Nature Genetics synteny list
#SYN=~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/Mo17_Table2_gene_list/1.B73_Syntenic_genes/19.B73-sytenic_gene.txt
SYN=~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/B73v4.36_Mo17CAU_synfind_synmap_SASoutput_avn.tsv

for GENO in B73 C123 Hp301 Mo17 NC338; do
    for TRT in Amb Ele; do
        echo "${GENO} ${TRT}"
        LOG=${OUTD}/${GENO}_${TRT}.log
        echo "${GENO} ${TRT} Log
" > ${LOG}

        ## Merge mean mapped read vlaues with 12604 and B73-Mo17 synteny
        python $PROJ/scripts/merge_b73_mo17_syntenic_cvr_cnt.py \
            -s ${SYN} \
            -g ${GENE} \
            -b ${INB}/cvrg_shrtRead_reads_in_region_${GENO}_sbys.csv \
            -m ${INM}/cvrg_shrtRead_reads_in_region_${GENO}_sbys.csv \
            -t ${TRT} \
            -v "sum_reads_in_region" \
            -1 \
            -p ${OUTD}/${GENO}_${TRT}_12604_cvr_cnt_1to1 \
            -c ${OUTD}/${GENO}_${TRT}_12604_1to1_reads_in_region_gene_counts.csv \
            >> ${LOG}

        ## Plot scatter plot of means for one-to-one pairs
        ## NOTE: the label names are written so that python will format the log10 with subscripts
        python ${SCRIPTS}/plot_scatter_from_CSV.py \
            -i ${OUTD}/${GENO}_${TRT}_12604_cvr_cnt_1to1_sum_reads_in_region.csv \
            -x sum_reads_in_region_${GENO}_${TRT}_mapB73 \
            -y sum_reads_in_region_${GENO}_${TRT}_mapMo17 \
            -xlab "\$Log_{10}\$ Sum Mapped Reads to B73" \
            -ylab "\$Log_{10}\$ Sum Mapped Reads to Mo17" \
            -r \
            -s log10 \
            -t png \
            -o ${OUTD}/${GENO}_${TRT}_12604_cvr_cnt_1to1_sum_reads_in_region.png

    done
done

## Combine gene counts files
SUMMARY=${OUTD}/all_12604_1to1_reads_in_region_gene_counts.csv
HEADER=$(head -1 ${OUTD}/B73_Amb_12604_1to1_reads_in_region_gene_counts.csv | \
        awk '{print "genotype,treatment,"$0}')

for GENO in B73 C123 Hp301 Mo17 NC338; do
    for TRT in Amb Ele; do
        awk -v geno=${GENO} -v trt=${TRT} 'NR!=1{print geno","trt","$0}' \
            ${OUTD}/${GENO}_${TRT}_12604_1to1_reads_in_region_gene_counts.csv
    done
done | awk -v header=${HEADER} 'BEGIN{print header}{print $0}' > ${SUMMARY}
rm ${OUTD}/*_*_12604_1to1_reads_in_region_gene_counts.csv
