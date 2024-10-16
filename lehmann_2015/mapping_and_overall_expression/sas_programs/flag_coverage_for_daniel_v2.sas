libname SEM '!MCLAB/sem_gof/sasdata';
libname SEMLOCAL '/home/jfear/storage/sem_gof/sasdata';
libname CEGS '!MCLAB/cegs_sergey/sas_data';

/* add on the incomplete data that I have droped from future analysis */

data tmp;
    set SEM.complete_ccfus_stacked SEMLOCAL.incomp_ccfus_stacked_all;
    run;

/* Using the data I imported for SEMs I need to create two flags */
    * flag_ge_one_read and flag_ge_hundred_read;
data flag;
    set tmp;
    if reads_in_region ge 1 then flag_ge_one_read = 1; else flag_ge_one_read = 0;
    if reads_in_region ge 100 then flag_ge_hundred_read = 1; else flag_ge_hundred_read = 0;
    genotype_rep = trim(line) || '_' || trim(mating_status) || trim(rep);
    keep line mating_status rep fusion_id reads_in_region flag_ge_one_read flag_ge_hundred_read genotype_rep;
    run;

proc sort data=flag;
    by genotype_rep;
    run;

proc means data =flag noprint;
    by genotype_rep;
    output out=sums sum(flag_ge_one_read)=num_fusion_ge_one_read sum(flag_ge_hundred_read)=num_fusion_ge_hundred_read /autoname;
    run;

data percents;
    retain genotype_rep _freq_ num_fusion_ge_one_read per_fusion_ge_one_read num_fusion_ge_hundred_read per_fusion_ge_hundred_read;
    set sums;
    rename _freq_ = total_num_fusions;
    per_fusion_ge_one_read = num_fusion_ge_one_read / _freq_ * 100;
    per_fusion_ge_hundred_read = num_fusion_ge_hundred_read / _freq_ * 100;
    drop _type_;
    run;

/* import Daniels table */
proc import datafile='/home/jfear/mclab/cegs_sergey/original_data/documentation/daniels_list_v2.csv' out=daniel dbms=csv replace;
    getnames=yes;
    run;


/* merge two datasets together */
proc sort data=daniel;
    by genotype_rep;
    run;
proc sort data=percents;
    by genotype_rep;
    run;

data merged oops1 oops2;
    merge daniel (in=in1) percents (in=in2);
    by genotype_rep;
    if in1 and in2 then output merged;
    if in1 and not in2 then output oops1;
    if in2 and not in1 then output oops2;
    run;

data oops1;
    set oops1;
    keep genotype_rep;
    run;

proc export data=oops1 outfile='/home/jfear/mclab/cegs_sergey/reports_external/samples_only_in_daniels_list.csv' dbms=csv replace;
    putnames =yes;
    run;

proc sort data=merged;
    by num_fusion_ge_hundred_read;
    run;

proc export data=merged outfile='/home/jfear/mclab/cegs_sergey/reports_external/coverage_count_flags_v2.csv' dbms=csv replace;
    putnames =yes;
    run;

proc sgplot data=merged;
    histogram num_fusion_ge_hundred_read;
    run;
