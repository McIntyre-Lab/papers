/*
 * REVISIONS: 12/21/2011
 *            - Removed Freq Results, now published in NOTEBOOK
 */

proc sort data=fru.all_coverage_counts_miss_w_key;
    by trt;
    run;

data fru.flag_no_trt_miss;
    set fru.all_coverage_counts_miss_w_key;
    if trt = "" then flag_no_trt = 1;
    else flag_no_trt = 0; 
    keep fusion_id sample_id flag_no_trt;
    run;

proc freq data=fru.flag_no_trt_miss;
    tables flag_no_trt;
    run;

data flag_dsx;
    set fru.all_coverage_counts_miss_w_key;
    where trt ? 'dsx'; 
    flag_dsx = 1;
    keep fusion_id sample_id flag_dsx;
    run;

proc sort data=fru.all_coverage_counts_miss_w_key;
    by fusion_id sample_id;
    run;

proc sort data=flag_dsx;
    by fusion_id sample_id;
    run;

data flag_dsx2;
    merge fru.all_coverage_counts_miss_w_key (in=in1) flag_dsx (in=in2);
    by fusion_id sample_id;
    run;

data fru.flag_dsx_miss;
    set flag_dsx2;
    if flag_dsx ne 1 then flag_dsx=0;
    keep fusion_id sample_id flag_dsx;
    run;




