### Supplementary cluster results ###

Descriptions of directory contents and associated documentation files are listed below
NOTE: All documentation can be found https://github.com/McIntyre-Lab/papers/tree/master/nanni_maize_2022/documentation


## Directory contents ##

* cluster_map_id_cov_counts/ :
    Tables containing the number and proportion of clusters mapped to
        B73 v4 (b73), Mo17 Yan (mo17_yan), or Mo17 Cau (mo17_cau) reference
        genomes with at least 95, 90, 85, or 80 percent sequence identity
        based on the alignment. Also included are the number and proportion
        of clusters with an alignment length that is at least 95, 90, 85, or
        80 percent of the total cluster sequence length.
        See documentation below for full column header descriptions.
    Documentation = map_ccs_reads.xlsx (see row 30)

* mapped_cluster_2_NAM_pangene_counts.txt :
    Counts comparing mapped clusters from each genotype-treatment sample to
        the NAM pan-gene classifications (Core Gene, Dispensable Gene,
        Near-Core Gene, Private Gene).
    Documentation = cluster_mapping_evaluation.xlsx (see row 14)

* compare_clusters_in_B73_Mo17_synteny/ :
    Counts for each of the Mo17 samples
        (*_b73_mo17CAU_cluster_compare_counts.txt) and all 3 Mo17 samples
        combined (combined_counts.txt) comparing clusters mapping to B73 or
        Mo17 CAU within B73-Mo17 syntenic gene pairs.
    Documentation = cluster_mapping_evaluation.xlsx (see row 20)

* ignored_clusterLength_distributions/ :
    PNG histogram plots of cluster sequence lengths ignored due to low
        coverage (<99%) (*_unCollapsed_low_cov_b73_ref_seq_length_hist.png)
        or low identity (<95%)
        (*_unCollapsed_low_identity_b73_ref_seq_length_hist.png) in mapping to
        the B73 v4 reference genome for each long read sample.
    Documentation = cluster_mapping_evaluation.xlsx (see row 32)

* cluster_fastqc_output/ :
    FASTQC output of cluster sequnces ignored due to low coverage (<99%)
        (*_unCollapsed_low_cov_b73_ref_fastqc.html) or ignored due to low
        identity (<95%) (*_unCollapsed_low_ident_b73_ref_fastqc.html) in
        mapping to the B73 v4 reference genome for each sample. Clusters
        passing the thresholds are collapsed within each sample to form
        assembled transcripts, also evaluated with FASTQC
        (*.collapsed.rep_fastqc.html).
    Documentation = cluster_mapping_evaluation.xlsx (see row 42)

* ignored_cluster_blast_results/ :
    Summary of BLAST results for the clusters idnored by Cupcake ToFU2 due
        to low identity (<95% of reference sequence identity) or low coverage
        (<99% of cluster sequence) when mapped to the B73 v4 reference genome
        compared to the NCBI non-redundant nucleotide database
        (summary_counts_2_b73.csv). See documentation below for full column
        header descriptions. Additional results include the C123, Hp301, and
        NC338 ignored clusters compared to the assembled transcripts of B73
        and Mo17 samples using BLAST (summary_counts_2_b73_vs_b73_mo17_clusters.csv).
        FASTA files used for input to BLAST, BLAST output files, and best hits
        for each cluster (as defined by the largest bitscore and smallest
        e-value) are available on zenodo ().
    Documentation = cluster_mapping_evaluation.xlsx (see row 32 and 38)

