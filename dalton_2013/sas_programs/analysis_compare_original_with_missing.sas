/* Compare the original analysis to the new analysis with the missing added 
 * Since we noticed that the Male FruM(A) rep 3 was truncated, we are
 * wanting to make sure that we check the fdr flags and see how many change
 * by including these additional (~800k reads). 
 */

proc sort data=FRU.flag_fdr_contrast_by_fusion;
by fusion_id;
run;

proc sort data=FRU.flag_fdr_contrast_by_fusion_miss;
by fusion_id;
run;

data merged oops;
    merge FRU.flag_fdr_contrast_by_fusion (in=in1) FRU.flag_fdr_contrast_by_fusion_miss (in=in2);
    by fusion_id;
    if in1 and in2 then output merged;
    else output oops;
    run;


data merged2 ;
    merge FRU.flag_fdr_contrast_by_fusion (in=in1) FRU.flag_fdr_contrast_by_fusion_miss (in=in2);
    by fusion_id;
    run;

proc contents data=merged2 varnum;
run;

proc freq data=merged2;
table flag_fdr_p_contrast_11_20*flag_fdr_p_miss_11_20;
table flag_fdr_p_contrast_12_20*flag_fdr_p_miss_12_20;
table flag_fdr_p_contrast_9_20*flag_fdr_p_miss_9_20;
run;

data tmp;
    set merged2;
    keep fusion_id flag_fdr_p_contrast_11_20 flag_fdr_p_miss_11_20;
    run;
