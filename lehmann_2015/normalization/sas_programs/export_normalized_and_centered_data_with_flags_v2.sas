/* Combine Raw, normalized, centered datasets */
    * All coverage Count data ;
    data raw;
        set CEGLOCAL.ccfus_stack;
        drop sample;
        run;

    * Noramlized coverage count data;
    data norm; 
        set CEGLOCAL.norm_basic_stats_m CEGLOCAL.norm_basic_stats_v;
        drop mapped_reads read_length region_length region_depth reads_in_region apn rpkm sample;
        run;

    * Centered and coverage counts;
    data centered;
        set CEGLOCAL.ccfus_norm_centered;
        drop apn sample uq_apn uq_ff log_uq_apn flag_raleigh;
        run;

    * Sort and Merge;
    proc sort data=raw;
        by fusion_id line mating_status rep;
        run;

    proc sort data=norm;
        by fusion_id line mating_status rep;
        run;

    proc sort data=centered;
        by fusion_id line mating_status rep;
        run;

    data alldata;
        merge raw (in=in1) norm (in=in2) centered (in=in3);
        by fusion_id line mating_status rep;
        run;

/* Merge On Various Flags */

    * flag_drop_fusion;
    proc sort data=CEGS.flag_drop_fusion;
        by fusion_id mating_status;
        run;

    proc sort data=alldata;
        by fusion_id mating_status;
        run;

    data alldata_flag1;
        merge alldata (in=in1) CEGS.flag_drop_fusion (in=in2);
        by fusion_id mating_status;
        run;

    * flag_raleigh;
    proc sort data=CEGS.flag_raleigh;
        by line;
        run;

    proc sort data=alldata_flag1;
        by line;
        run;

    data alldata_flag2;
        merge alldata_flag1 (in=in1) CEGS.flag_raleigh (in=in2);
        by line;
        run;

    * flag_mahalanobis_outlier, flag_sample_on, flag_dataset_of_origin;
    proc sort data=CEGS.flag_mahalanobis_outlier;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_sample_on;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_dataset_of_origin;
        by line mating_status rep;
        run;

    proc sort data=CEGS.detected_exon_cnts;
        by line mating_status rep;
        run;

    proc sort data=alldata_flag2;
        by line mating_status rep;
        run;

    data alldata_flag3;
        merge alldata_flag2 (in=in1) CEGS.flag_mahalanobis_outlier (in=in2) CEGS.flag_sample_on (in=in3) CEGS.flag_dataset_of_origin (in=in4) CEGS.detected_exon_cnts (in=in5);
        by line mating_status rep;
        run;

/* Merge on gene information */
    data dmel;
        set DMEL.fb551_si_fusions_unique_flagged;
        keep fusion_id symbol_cat fbgn_cat constitutive common alternative;
        run;
        
    proc sort data=dmel;
        by fusion_id;
        run;

    proc sort data=alldata_flag3;
        by fusion_id;
        run;

    data merged;
        merge alldata_flag3 (in=in1) dmel (in=in2);
        by fusion_id;
        if in1;
        run;

/* Re-organize columns for export */
    proc contents data=merged varnum;
    run;

    data merged2;
        retain line mating_status rep flag_raleigh flag_incomplete flag_partial
        flag_complete fusion_id mapped_reads read_length region_length region_depth
        reads_in_region apn rpkm detected_exon_cnts flag_sample_on flag_drop_fusion
        flag_mahalanobis_outlier sum_mapped q3 median junc_mapped_reads total_mapped_reads uq_apn
        uq_ff log_uq_apn mean_log_uq_apn median_log_uq_apn uq_log_uq_apn
        mean_log_uq_center median_log_uq_center uq_log_uq_center symbol_cat FBgn_cat
        Constitutive Common Alternative;
        set merged;
        run;

/* Export */
    proc export data=merged2 outfile='!MCLAB/cegs_sergey/pipeline_output/rawAPN_uqNormCenter_plusFlags_20140518.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Clean UP */
    proc datasets nolist;
        delete alldata;
        delete alldata_flag1;
        delete alldata_flag2;
        delete alldata_flag3;
        delete dmel;
        delete centered;
        delete merged;
        delete merged2;
        delete norm;
        delete raw;
        run; quit;
