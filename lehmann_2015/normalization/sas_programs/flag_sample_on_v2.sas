/********************************************************************************
* First determine which samples actually have enough expression. I will say a
* sample is on if it has at least 2k exonic regions that have an APN >= 5.
********************************************************************************/

libname OE '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_counts/sas_data';
libname DESIGN '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/sas_data';
libname NORM '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';

data alldata;
    set OE.ccfus_stack;
    if APN > 0 then flag_expressed = 1; else flag_expressed = 0;
    if APN ge 5 then flag_expressed5 = 1; else flag_expressed5 = 0;
    run;

proc sort data=alldata;
    by line mating_status rep;
    run;

proc means data=alldata noprint;
    by line mating_status rep;
    output out=sums sum(flag_expressed)=cnt_apn_gt_0 sum(flag_expressed5)=cnt_apn_gt_5;
    run;

/* Look at the numbers */
proc sort data=sums;
    by cnt_apn_gt_0;
    run;

/* Now lets combine and make lots of flags */
proc sort data=sums;
    by line mating_status rep;
    run;

proc sort data=DESIGN.flag_dataset_of_origin;
    by line mating_status rep;
    run;

data merged;
    merge sums (in=in1) DESIGN.flag_dataset_of_origin (in=in2);
    by line mating_status rep;
    if cnt_apn_gt_0 >= 29300 then flag_apn_gt_0_gt_29k = 1; else flag_apn_gt_0_gt_29k = 0;
    if flag_complete = 1 then flag_data = 'complete';
    else if flag_partial = 1 then flag_data = 'partial';
    else flag_data = 'incomplete';
    drop _type_ ;
    run;

proc export data=merged outfile="/tmp/sums.csv" dbms=csv replace;
    putnames=yes;
    run;

/* Cound not get colors to work in proc iml */
data _null_;
    call system('Rscript $MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/r_programs/plot_count_cutoff.R /tmp/sums.csv');
    run;

data NORM.flag_sample_on;
    set merged;
    if flag_apn_gt_0_gt_29k=1 then flag_sample_on = 1; 
    else flag_sample_on = 0;
    keep line mating_status rep flag_sample_on;
    run;

proc freq data=NORM.flag_sample_on;
    table flag_sample_on;
    run; *dropped 97 samples;

proc datasets nolist;
    delete alldata sums merged ;
    run; quit;
