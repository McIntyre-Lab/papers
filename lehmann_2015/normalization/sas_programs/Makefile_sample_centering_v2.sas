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
    * INPUT: CEGLOCAL.ccfus_norm_stack_m 
    *        CEGLOCAL.ccfus_norm_stack_v
    *
    * DATASET: CEGS.ccfus_norm_centered
    *          CEGLOCAL.ccfus_norm_centered (For Faster Access)
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/normalize_center_by_sample.sas';

/* Distribution Plots */
    * Plot the distrubtion of the Raw APN, log_uq_apn, and centered values
    * 
    * INPUT: CEGLOCAL.ccfus_norm_centered
    * 
    * RSCRIPT: $MCLAB/cegs_sergey/r_programs/plot_normalization_distributions.R
    *
    * FIGURES: !MCLAB/cegs_sergey/reports/line_normalization/Mated_uq_log_uq_centered_boxplot_v2.png
    *          !MCLAB/cegs_sergey/reports/line_normalization/Virgin_uq_log_uq_centered_boxplot_v2.png
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/plot_normalization_distributions.sas';

/* Create Side-by-side dataset for centered data */
    * Want to look at some correlations and whatnot, so I need a side-by-side
    * dataset.
    *
    * INPUT: CEGLOCAL.ccfus_norm_centered
    *
    * DATASET: CEGLOCAL.norm_sbs
    *          CEGLOCAL.center_sbs
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/centered_stack_to_SBS.sas';

/* Use JMP to calculate Mahalanobis distance */
    * I used JMP to calculate Mahalanobis distance and output a datasets. 
    *
    * Genomics > Expression > Quality Control > Correlation and Principal Component Analysis
    *
    * INPUT: CEGLOCAL.center_sbs
    * 
    * DATASET: CEGS.mahalanobis_distance
    ;

/* Flag Outliers */
    * Create flags of the outliers (mahalanobis distance >= 4) and create
    * boxplots. If an entire line is an outlier then drop.
    *
    * All flaged as outlier 
    *
    *     r486 M
    *     r491 M
    *     w23 M
    *     w67 M
    *     w82 M
    *     w82 V
    *
    * INPUT: CEGS.mahalanobis_distance
    * 
    * RSCRIPT: $MCLAB/cegs_sergey/r_programs/plot_mahalanobis_outliers.R 
    *
    * DATASET: CEGS.flag_mahalanobis_outlier 
    * 
    * FILES: $MCLAB/cegs_sergey/reports/mahalanobis_outlier.png
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_outlier.sas';
