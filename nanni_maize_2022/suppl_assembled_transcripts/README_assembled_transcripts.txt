### Supplementary assembled transcript results ###

Descriptions of directory contents and associated documentation files are listed below
NOTE: All documentation can be found https://github.com/McIntyre-Lab/papers/tree/master/nanni_maize_2022/documentation


## Directory contents ##

* assembXcrpt_length_distributions/ :
    PNG histogram plots of assembled transcript sequence lengths following
        TOFU2 Cupcake collapse of cluster sequences passing the thresholds
        of coverage (>=99%) and identity (>=95%) in mapping to the B73 v4
        reference genome for each long read sample
        (*_pass_tofu_filter_b73_ref_seq_length_hist.png).
    Documentation = cluster_mapping_evaluation.xlsx (see row 40)

* DE_GO_results/ :
    For each Gene Ontology (GO) category (biologcal process, cellular component,
        molecular function), combined GO enrichment results from 4 different
        enrichment tests of DE groups (comb_GO_gse_tappas_DE_*_sig_terms.csv) from:
        1) RSEM quantification of B73 v4 reference genes with assembled transcripts (rsemB73map),
        2) genome alignments of all B73 v4 reference genes (ccB73map),
        3) genome alignments of B73 v4 reference genes within 1-to-1 B73-Mo17 syntenic pairs (ccB73mapBMo1to1),
        and 4) genome alignments of Mo17 CAU reference genes within 1-to-1 B73-Mo17 syntenic pairs (ccMo17mapBMo1to1).
        Also included for each GO category is a count of any shared GO terms aross
        the 4 tests (*_identical_all_4_enrichments.csv), pairwise comparison of enrichment
        tests rsemB73map vs. ccB73map (*_identical_rsemB73map_ccB73map_enrichments.csv),
        pairwise comparison of enrichment tests ccB73map vs. ccB73mapBMo1to1
        (*_identical_ccB73map_ccB73mapBMo1to1_enrichments.csv), and pairwise comparison of
        enrichment tests ccB73mapBMo1to1 and ccMo17mapBMo1to1 (*__identical_ccB73mapBMo1to1_ccMo17mapBMo1to1_enrichments.csv).
        See documentation below for full column header descriptions.
    Documentation = shrt_read_alt_analysis.xlsx (see row 50)

* B73_vs_Mo17_CAU_mapped_read_scatter/ :
    For each genotype-treatment, PNG files of scatter plots comparing the number of
        mapped reads in exonic regions of genes (sum of reads in region) when mapped
        to B73 v4 reference genome or Mo17 CAU reference genome on a log10 scale
        for each 1-to-1 B73-Mo17 syntenic gene pair. Pearson's correlation coefficient (r)
        is indicated for each comparion.
    Documentation = shrt_read_alt_analysis.xlsx (see row 54)
