/********************************************************************************
* This script combines our result table with fusion level flags for if a fusion
* was considred induced or repressed according to our criteria.
********************************************************************************/

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';

%macro sort_by_fusion(indata);
    proc sort data=FRU.&indata;
        by fusion_id;
        run;
%mend;

%sort_by_fusion(results_by_fusion);
%sort_by_fusion(frua_male_fusions_ind_rep);
%sort_by_fusion(frua_female_fusions_ind_rep);
%sort_by_fusion(frub_male_fusions_ind_rep);
%sort_by_fusion(frub_female_fusions_ind_rep);
%sort_by_fusion(fruc_male_fusions_ind_rep);
%sort_by_fusion(fruc_female_fusions_ind_rep);
%sort_by_fusion(null_male_fusions_ind_rep);

data FRU.results_with_flag_ind_rep;
    retain fusion_id symbol_cat flag_fail_normality 
    flag_fruA_male_ind flag_fruA_male_rep flag_fruA_female_ind
    flag_fruA_female_rep flag_fruB_male_ind flag_fruB_male_rep
    flag_fruB_female_ind flag_fruB_female_rep flag_fruC_male_ind
    flag_fruC_male_rep flag_fruC_female_ind flag_fruC_female_rep 
    flag_null_male_ind flag_null_male_rep;
    merge FRU.results_by_fusion FRU.frua_male_fusions_ind_rep
    FRU.frua_female_fusions_ind_rep FRU.frub_male_fusions_ind_rep
    FRU.frub_female_fusions_ind_rep FRU.fruc_male_fusions_ind_rep
    FRU.fruc_female_fusions_ind_rep FRU.null_male_fusions_ind_rep;
    by fusion_id;
    drop Genes_per_fusion Exons_per_fusion FBgns_per_fusion FBpps_per_fusion FBtrs_per_fusion CHROM
    START END Exon_Gene_ID_cat Sequence_loc_cat exon_ID_cat Exon_Name_cat FBtrs_per_exon_cat FBpp_cat FBtr_cat
    max_fbtr_per_gene_symbol min_fbtr_per_gene_symbol_cat max_fbtr_per_gene_symbol_cat;
    run;
