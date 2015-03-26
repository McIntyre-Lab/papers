/*
 * This script will go ahead and drop fusion_id*sample_id from the main dataset
 */

proc sort data=fru.all_coverage_counts_miss_w_key;
    by fusion_id sample_id;
    run;

proc sort data=fru.flag_no_trt_miss;
    by fusion_id sample_id;
    run;

proc sort data=fru.flag_dsx_miss;
    by fusion_id sample_id;
    run;

data short_list oops;
    merge fru.all_coverage_counts_miss_w_key (in=in1) fru.flag_no_trt_miss (in=in2) fru.flag_dsx_miss (in=in3);
    by fusion_id sample_id;
    if in1 and in2 and in3 then output short_list;
    else output oops;
    run;

data short_list2;
    set short_list;
    where flag_no_trt = 0 and flag_dsx = 0;
    run;

