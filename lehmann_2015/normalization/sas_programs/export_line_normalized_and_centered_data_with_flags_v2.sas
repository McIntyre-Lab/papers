/* Combine Raw, normalized, centered datasets */
    * Noramlized coverage count data;
    data norm; 
        set CEGLOCAL.line_norm_basic_stats_m CEGLOCAL.line_norm_basic_stats_v;
        *drop sum_mapped_reads read_length region_length region_depth reads_in_region apn uq_ff q3 median;
        run;

    * Centered and coverage counts;
    data centered;
        set CEGLOCAL.line_norm_centered;
        drop apn uq_apn uq_ff log_uq_apn mean_log_uq_apn median_log_uq_apn uq_log_uq_apn;
        run;

    * Sort and Merge;
    proc sort data=norm;
        by fusion_id line mating_status ;
        run;

    proc sort data=centered;
        by fusion_id line mating_status ;
        run;

    data alldata;
        merge norm (in=in1) centered (in=in2);
        by fusion_id line mating_status ;
        run;

/* Merge On Various Flags */

    * flag_raleigh;
    proc sort data=CEGS.flag_raleigh;
        by line;
        run;

    proc sort data=alldata;
        by line;
        run;

    data alldata_flag1;
        merge alldata (in=in1) CEGS.flag_raleigh (in=in2);
        by line;
        if in1;
        run;

/* Merge on gene information */
    data dmel;
        set DMEL.fb551_si_fusions_unique_flagged;
        keep fusion_id symbol_cat fbgn_cat constitutive common alternative;
        run;
        
    proc sort data=dmel;
        by fusion_id;
        run;

    proc sort data=alldata_flag1;
        by fusion_id;
        run;

    data merged;
        merge alldata_flag1 (in=in1) dmel (in=in2);
        by fusion_id;
        if in1;
        run;

/* Re-organize columns for export */
    proc contents data=merged varnum;
    run;

    data merged2;
        retain line mating_status fusion_id Constitutive Common Alternative
        sum_region_depth sum_mapped uq_apn log_uq_apn mean_log_uq_center
        median_log_uq_center uq_log_uq_center flag_raleigh symbol_cat FBgn_cat;
        set merged;
        run;

/* Export */
    *proc export data=merged2 outfile='!MCLAB/cegs_sergey/pipeline_output/line_level_uqNormCenter_plusFlags_20140518.csv' dbms=csv replace;
    proc export data=merged2 outfile='!MCLAB/cegs_sergey/pipeline_output/line_level_uqNormCenter_plusFlags_w_uqff_for_david.csv' dbms=csv replace;
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
