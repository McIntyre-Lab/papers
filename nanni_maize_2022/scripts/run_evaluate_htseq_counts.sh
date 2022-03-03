#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/sqanti_post_filter_b73
OUTD=$PROJ/htseq_output_counts
    mkdir -p ${OUTD}
ROZ=${OUTD}/roz_eval_htseq
    mkdir -p ${ROZ}
LOG=${OUTD}/evaluate_htseq_counts.log


if [[ ! -e ${OUTD}/combined_htseq_count.tsv ]]; then
    scp adalena.nanni@hpg.rc.ufl.edu:/blue/mcintyre/share/maize_ainsworth/htseq_quantify_genes_loci/htseq_output_counts/combined_htseq_count.tsv ${OUTD}/.
fi

## Average gene/loci expression across samples in each genotype/treatment group
## All genes/loc
echo "All genes = $(awk 'NR!=1' ${OUTD}/combined_htseq_count.tsv)" > ${LOG}
python $PROJ/scripts/flag_on_off_rsem_02avn.py \
    -i ${OUTD}/combined_htseq_count.tsv \
    -m 0 \
    -t gene \
    -u reads \
    -g 0 -g 3 \
    -p 0.05 \
    --plot \
    -o ${OUTD}/flag_on_off_htseq \
    >> ${LOG}

## Novel/fusion/antisense loci only
awk -F "\t" 'NR==1 || $1~"_"' ${OUTD}/combined_htseq_count.tsv \
    > ${ROZ}/novel_loci.tsv
echo "$(awk 'NR!=1' ${ROZ}/novel_loci.tsv | wc -l) novel/fusion/antisense loci" >> ${LOG}
echo -e "\t$(awk -F "\t" 'NR!=1 && $1~"AS"' ${OUTD}/combined_htseq_count.tsv | wc -l) antisense loci" >> ${LOG}
echo -e "\t$(awk -F "\t" 'NR!=1 && $1!~"novel" && $1~"_"' ${OUTD}/combined_htseq_count.tsv | wc -l) fusion loci" >> ${LOG}
echo -e "\t$(awk -F "\t" 'NR!=1 && $1!~"AS" && $1~"novel"' ${OUTD}/combined_htseq_count.tsv | wc -l) other novel loci" >> ${LOG}
python $PROJ/scripts/flag_on_off_rsem_02avn.py \
    -i ${ROZ}/novel_loci.tsv \
    -m 0 \
    -t gene \
    -u reads \
    -g 0 -g 3 \
    -p 0.05 \
    --plot \
    -o ${OUTD}/flag_on_off_htseq_novel_loci \
    >> ${LOG}

## Evaluate htseq counts for novel loci
python $PROJ/scripts/evaluate_htseq_counts.py \
    -f ${OUTD}/flag_on_off_htseq_novel_loci_0_full_table.tsv \
    -p ${OUTD}/novel_loci_htseq \
    >> ${LOG}

## Make fasta file of PB transcripts
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
gffread -w ${ROZ}/sqanti_post_filter_b73.fasta \
    -g ${REF} \
    $PROJ/sqanti_post_filter_b73/sqanti_b73_filtered_corrected_associated_gene.gtf
## Convert from multifasta to single line fasta
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
    ${ROZ}/sqanti_post_filter_b73.fasta \
    > ${ROZ}/sqanti_post_filter_b73.single.fasta

## Get PBid values for the 3 potentially B73 exclusive loci with > 5 reads in both treatments
ALLID=""
for LOCI in novelGene_372 novelGene_344 novelGene_769; do
    PBID=$(awk -F "\t" -v loci=${LOCI} '$7==loci{total=total"\t"$1}END{print total}' $PROJ/sqanti_post_filter_b73/SQANTI_classification.txt)
    ALLID=$(echo -e "${ALLID}\t${PBID}")
    echo -e "novelGene from SQANTI: ${LOCI}\n\tPBid values:${PBID}"
done

head -1 $PROJ/merge_samples_b73/all_samples.chained_ids.txt \
    > ${OUTD}/chained_sample_subset.tsv
touch ${OUTD}/novel_loci_htseq_B73_only_novel_loci_ge5reads.fasta
for ID in ${ALLID}; do
    awk -v id=${ID} -F "\t" '$1==id' \
        $PROJ/merge_samples_b73/all_samples.chained_ids.txt \
        >> ${OUTD}/chained_sample_subset.tsv
    grep -A 1 ">${ID}" ${ROZ}/sqanti_post_filter_b73.single.fasta \
        >> ${OUTD}/novel_loci_htseq_B73_only_novel_loci_ge5reads.fasta
done

rm -r ${ROZ}

## Copy fasta to share
#scp ${OUTD}/novel_loci_htseq_B73_only_novel_loci_ge5reads.fasta \
#    adalena.nanni@hpg.rc.ufl.edu:/blue/mcintyre/share/maize_ainsworth/htseq_quantify_genes_loci/.
