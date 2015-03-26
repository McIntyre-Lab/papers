/*
 * This script will go ahead and drop fusion_id*sample_id from the main dataset
 */

libname fru '!MCLAB/Fru_network/sasdata';

proc sort data=fru.all_coverage_counts_with_key;
    by fusion_id sample_id;
    run;

proc sort data=fru.flag_no_trt;
    by fusion_id sample_id;
    run;

proc sort data=fru.flag_dsx;
    by fusion_id sample_id;
    run;

data short_list oops;
    merge fru.all_coverage_counts_with_key (in=in1) fru.flag_no_trt (in=in2) fru.flag_dsx (in=in3);
    by fusion_id sample_id;
    if in1 and in2 and in3 then output short_list;
    else output oops;
    run;

data short_list2;
    set short_list;
    where flag_no_trt = 0 and flag_dsx = 0;
    run;

