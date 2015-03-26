/*
 * REVISIONS: 12/21/2011
 *            - Changed libname
 */
libname fru '!MCLAB/Fru_network/sasdata';

proc sort data = short_list2;
    by trt fusion_id;
    run;

/* calculate the variance for each trt*fusion */
proc means data=short_list2 noprint ;
    class trt fusion_id;
    var rpkm;
    output out=counts_vars var=var;
    run;

/* create a flag for each trt*fusion */
data flag_fusion_trt_no_var;
    set counts_vars;
    where _type_=3;
    if var=0 then flag_no_var =1;
    else flag_no_var = 0;
    run;

proc sort data=flag_fusion_trt_no_var;
    by fusion_id;
    run;

proc means data=flag_fusion_trt_no_var noprint;
    by fusion_id;
    output out=sum_flags_no_var sum(flag_no_var)=sum_flag_no_var;
    run;

data fru.flag_no_var;
    set sum_flags_no_var;
    if sum_flag_no_var > 0 then flag_no_var=1;
    else flag_no_var=0;
    keep fusion_id flag_no_var;
    run; *60291 obs;

proc freq data=fru.flag_no_var;
    tables flag_no_var;
    run;

