/********************************************************************************
* We are interested in the overlap between genes that were induced by fru
* overexpression and ko. I will do all pairwise compairsons between the
* different overexpression {A,B,C} and null.
********************************************************************************/

libname fru '!MCLAB/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Create ind rep no_DEG flags */
    * In order to do the comparison that we are wanting to do, I need to create new
    * flags (-1,0,1) where (-1) is repressed, (0) is no change, (1) is induced.
    *
    * INPUT: FRU.flag_x_induced_repressed_male
    *        FRU.flag_x_ind_rep_null_male
    *
    * DATASET: FRU.flag_male_status
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_enrich_flag_ind_rep_for_comparison.sas';

/* Fisher's Exact Test */
    * Calculate Fisher's exact test to see if there is an association.
    *
    * INPUT: FRU.flag_male_status
    *
    * OUTPUT: "!MCLAB/arbeitman/arbeitman_fru_network/reports_external/male_ind_rep_enrichment_&type..csv"
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_enrich_male_enrichments_tests_over_vs_null.sas';
