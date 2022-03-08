### Supplementary CCS read results ###

Descriptions of directory contents and associated documentation files are listed below
NOTE: All documentation can be found https://github.com/McIntyre-Lab/papers/tree/master/nanni_maize_2022/documentation


## Directory contents ##

* ccs_readLength_distributions/ :
    PNG histogram plots of CCS read lengths for each genotype (${GENOTYPE}_ccs_seq_length_hist.png).
        A histogram of B73 v4 reference transcript lengths is included for comparison (ref_b73_seq_length_hist.png).
    Documentation = map_ccs_reads.xlsx (see row 16)

* ccs_map_id_cov_counts/ :
    Tables containing the number and proportion of CCS reads mapped to
        B73 v4 (b73) Mo17 Yan (mo17_yan) or Mo17 Cau (mo17_cau) reference genomes
        with at least 95, 90, 85, or 80 percent sequence identity based on the alignment.
        Also included are the number and proportion of CCS reads with an alignment length
        that is at least 95, 90, 85, or 80 percent of the total CCS read length.
        See documentation below for full column header descriptions.
    Documentation = map_ccs_reads.xlsx (see row 30)

* ccs_FSM_ISM_readLength_plots/ :
    PNG histogram plots of read lengths for CCS reads that map to B73 v4 reference genome and are
        full-splice match (FSM) to a reference transcript (${GENO}_FSM_ccs_seq_length_hist.png). 
        A scatter plot of FSM and incomplete-splice match (ISM) CCS read lengths vs.
        the associated B73 v4 reference transcript length is also included (${GENO}_FSM_ISM_vs_ref_length.png).
    Documentation = map_ccs_reads.xlsx (see row 40)

