/********************************************************************************
* First determine which samples actually have enough expression. I will say a
* sample is on if it has at least 2k exonic regions that have an APN >= 5.
********************************************************************************/

data alldata;
    set CEGLOCAL.platek_ccfus_stack CEGLOCAL.ccfus_stack ;
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

/* Now lets make some flags */
data CEGS.flag_sample_onk;
    set sums;
    if cnt_apn_gt_0 >= 29300  then flag_sample_on = 1; 
    else flag_sample_on = 0;
    keep line mating_status rep flag_sample_on;
    run;

proc freq data=CEGS.flag_sample_onk;
    table flag_sample_on;
    run; *dropped 104 samples;

proc datasets nolist;
    delete alldata sums merged ;
    run; quit;
