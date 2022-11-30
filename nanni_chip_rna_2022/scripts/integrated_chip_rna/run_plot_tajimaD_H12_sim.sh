#!/bin/sh

export PATH=$HOME/conda/.envs/mylab/bin:$PATH
#export PATH=/Users/adalena/opt/anaconda3/envs/mylab/bin:$PATH

## Make plots of distributions for Tajima's D and H12 for D. simulans

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
TD=~/mclab/SHARE/McIntyre_Lab/useful_dsim_data/gene_lists/tajimaD_h12_signor_2018/filter_nolab_tajd_1kb.Tajima.D
H12=~/mclab/SHARE/McIntyre_Lab/useful_dsim_data/gene_lists/tajimaD_h12_signor_2018
ANNOT=~/mclab/SHARE/McIntyre_Lab/useful_dsim_data/flybase202/dsim_annotation/fbgn2coord.csv
ORTHO=$PROJ/supp_files/dmel_dsim_ortholog_chip_rna_flags.csv

OUTD=$PROJ/figures/tajima_D
   mkdir -p ${OUTD}
ROZ=${OUTD}/roz_sim
   mkdir -p ${ROZ}

## Make BED file of gene coordinates
awk -F "," 'NR!=1{print $2"\t"$3"\t"$4"\t"$1}' \
    ${ANNOT} > ${ROZ}/sim_gene.bed

## Make BED file of Tajima's D windows from Signor, et al. 2018
awk 'NR!=1{print "Scf_"$1"\t"$2"\t"$2+1000"\t"$4}' \
    ${TD} > ${ROZ}/sim_tajimaD.bed

## Add header to intersect file
echo "gene_chr,gene_start,gene_end,gene_id,window_chr,window_start,window_end,tajimaD" \
    > ${OUTD}/dsim202_gene_tajimaD_intersect.csv

## Use bedtools intersect to associate Tajima's D values to genes
bedtools intersect \
    -a ${ROZ}/sim_gene.bed \
    -b ${ROZ}/sim_tajimaD.bed \
    -wa \
    -wb \
    | sed s/'\t'/','/g \
    >> ${OUTD}/dsim202_gene_tajimaD_intersect.csv

## Plot Tajima D boxplots
python $PROJ/scripts/integrated_chip_rna/plot_gene_group_distributions.py \
    -i ${OUTD}/dsim202_gene_tajimaD_intersect.csv \
    -v tajimaD \
    -n "Tajima's D" \
    -g ${ORTHO} \
    -p ${OUTD}/dsim202


## Make BED file of H12 peaks from Signor, et al. 2018
if [[ -e ${ROZ}/sim_h12.bed ]]; then
    rm ${ROZ}/sim_h12.bed
fi
for FILE in $(ls ${H12}/*_output.txt); do
    CHR=$(basename ${FILE} | awk -F "_" '{print $1}' | sed s/'chr'/'Scf_'/g)
    awk -v chr=${CHR} '{print chr"\t"$2"\t"$3"\t"$9"\t"$10"\t"$7}' \
        ${FILE} >> ${ROZ}/sim_h12.bed
done

## Add header to intersect file
echo "gene_chr,gene_start,gene_end,gene_id,peak_chr,peak_start,peak_end,H12,H2_over_H1,H1" \
    > ${OUTD}/dsim202_gene_H12_intersect.csv

## Use bedtools intersect to associate H12 values to genes
bedtools intersect \
    -a ${ROZ}/sim_gene.bed \
    -b ${ROZ}/sim_H12.bed \
    -wa \
    -wb \
    | sed s/'\t'/','/g \
    >> ${OUTD}/dsim202_gene_H12_intersect.csv

## Plot H12
python $PROJ/scripts/integrated_chip_rna/plot_gene_group_distributions.py \
    -i ${OUTD}/dsim202_gene_H12_intersect.csv \
    -v H12 \
    -n "H12" \
    -g ${ORTHO} \
    -p ${OUTD}/dsim202

## Plot H2/H1
python $PROJ/scripts/integrated_chip_rna/plot_gene_group_distributions.py \
    -i ${OUTD}/dsim202_gene_H12_intersect.csv \
    -v H2_over_H1 \
    -n "H2/H1" \
    -g ${ORTHO} \
    -p ${OUTD}/dsim202

## Plot H1
python $PROJ/scripts/integrated_chip_rna/plot_gene_group_distributions.py \
    -i ${OUTD}/dsim202_gene_H12_intersect.csv \
    -v H1 \
    -n "H1" \
    -g ${ORTHO} \
    -p ${OUTD}/dsim202

#rm -r ${ROZ}
