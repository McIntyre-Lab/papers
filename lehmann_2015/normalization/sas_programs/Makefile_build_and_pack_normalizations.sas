libname cegs '!MCLAB/cegs_sergey/sas_data';
libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Create Full Result Table with All Flags */
    * I want to create a full results table that includes all of the flags that
    * I have created so far as well as all of the normalized data.
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *        CEGLOCAL.ccfus_norm_basic_stats_m
    *        CEGLOCAL.ccfus_norm_basic_stats_v
    *        CEGLOCAL.ccfus_norm_centered
    *        CEGS.flag_fusion_on
    *        CEGS.flag_drop_fusion
    *        CEGS.detected_exon_cnts
    *        CEGS.flag_sample_on
    *        CEGS.flag_raleigh
    *        CEGS.flag_funky_boxplot
    *        CEGS.flag_dataset_of_origin
    *        DMEL.fb551_si_fusions_unique_flagged
    *
    * OUTFILE: !MCLAB/cegs_sergey/pipeline_output/rawAPN_uqNormCenter_plusFlags.csv
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/export_normalized_and_centered_data_with_flags_v2.sas';

