libname cegs '!MCLAB/cegs_sergey/sas_data';
libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Create Full Result Table with All Flags */
    * I want to create a full results table that includes all of the flags that
    * I have created so far as well as all of the normalized data.
    *
    * INPUT: CEGLOCAL.line_norm_basic_stats_m
    *        CEGLOCAL.line_norm_basic_stats_v
    *        CEGLOCAL.line_norm_centered
    *        CEGS.flag_raleigh
    *        CEGS.flag_dataset_of_origin
    *        DMEL.fb551_si_fusions_unique_flagged
    *
    * OUTFILE: !MCLAB/cegs_sergey/pipeline_output/line_level_rawAPN_uqNormCenter_plusFlags.csv
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/export_line_normalized_and_centered_data_with_flags_v2.sas';
