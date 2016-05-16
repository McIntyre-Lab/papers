/*
libname cegs '!MCLAB/cegs_sergey/sas_data';
libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/


/* Combine coverage files */
data coverage;
    set CEGS.ccfus_stack;
    run;

/* Summ across replicates */
proc sort data=coverage;
    by line mating_status fusion_id;
    run;

proc means data=coverage noprint;
    by line mating_status fusion_id;
    output out=sums sum(region_depth)=sum_region_depth;
    run;
    
data len;
    set coverage;
    keep fusion_id region_length;
    run;

proc sort data=len nodupkey;
    by fusion_id;
    run;

proc sort data=sums ;
    by fusion_id;
    run;

data merged1;
    merge sums len;
    by fusion_id;
    apn = sum_region_depth / region_length;
    drop _type_ _freq_;
    run;

/* Sort and Merge APN to bayesian machine */
proc sort data=merged1;
    by fusion_id line mating_status;
    run;

proc sort data=CEGS.flag_empirical_bayes;
    by fusion_id line mating_status;
    run;

data merged2 oops1 oops2;
    merge CEGS.flag_empirical_bayes (in=in1) merged1 (in=in2);
    by fusion_id line mating_status;
    if in1 and in2 then output merged2;
    if in1 and not in2 then output oops1;
    if in2 and not in1 then output oops2;
    run;
    
data forbox;
    set merged2;
    if flag_q4_AI =1 and flag_q5_AI =1 and flag_q6_AI= 1 then flag_all_AI = 1;
    else flag_all_AI = 0;
    keep line mating_status fusion_id flag_q4_AI flag_q5_AI flag_q6_AI flag_all_AI apn;
    run;

    proc export data=forbox outfile='!MCLAB/cegs_sergey/bayesian_analysis/empirical_bayesian_AI_vs_apn.csv' dbms=csv replace;
    putnames=yes;
    run;

