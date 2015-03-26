/*
 * REVISIONS: 12/21/2011
 *            - changed libname
 *            - [Line 14] changed logrpmk to rpkm
 *            - added freq statement
 *            12/23/2011
 *            - added a check step at end, both no var and zero mean are the same
 */

proc sort data=short_list2;
    by trt fusion_id;
    run; 

proc means data = short_list2 noprint;
    class trt fusion_id;
    var rpkm;
    output out=counts_means mean=mean;
    run;

/* create a dataset with the means for each trt variable sybs */
data fru.trt_means_miss;
    set counts_means;
    where _type_ = 3;
    if mean = 0 then flag_zero_trt_mean = 1;
    else flag_zero_trt_mean =0;
    run;

proc sort data=fru.trt_means_miss;
    by fusion_id;
    run;

proc means data=fru.trt_means_miss noprint;
    by fusion_id;
    output out=sum_flag_zero_mean sum(flag_zero_trt_mean)=sum_flag_zero_trt_mean;
    run;

data fru.flag_zero_mean_miss;
    set sum_flag_zero_mean;
    if sum_flag_zero_trt_mean > 0 then flag_zero_mean = 1;
    else flag_zero_mean =0;
    keep fusion_id flag_zero_mean;
    run;

proc freq data=fru.flag_zero_mean_miss;
    tables flag_zero_mean;
    run;

/* Checks 
data mean_check;
    set counts_means;
    where _type_ = 3 and mean = 0;
    run;
data var_check;
    set counts_vars;
    where _type_ = 3 and var = 0;
    run;
*/

