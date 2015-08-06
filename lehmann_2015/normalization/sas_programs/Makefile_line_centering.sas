libname cegs '!MCLAB/cegs_sergey/sas_data';
libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Centering */
    * From the log_uq_apn plots generated below, there is still a lot of
    * variation (coverage effects) that are not corrected for my the
    * uq-normalization. Now I am going to do a sample level centering. We are
    * trying 3 different centering strategies {Mean, Median, UQ}.
    *
    * INPUT: CEGLOCAL.line_norm_basic_stats_m 
    *        CEGLOCAL.line_norm_basic_stats_v
    *
    * DATASET: CEGS.line_norm_centered
    *          CEGLOCAL.line_norm_centered (For Faster Access)
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/normalize_center_by_line.sas';

/* Distribution Plots */
    * Plot the distrubtion of the centered values
    * 
    * INPUT: CEGLOCAL.ccfus_norm_centered
    * 
    * RSCRIPT: $MCLAB/cegs_sergey/r_programs/plot_normalization_distributions_line.R
    *
    * FIGURES: !MCLAB/cegs_sergey/reports/line_normalization/Mated_Virgin_uq_log_uq_centered_boxplot_line.pdf
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/plot_line_normalization_distributions.sas';
