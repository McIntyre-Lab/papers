/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

/* Get group means (APN, log_apn) */

proc sort data=mysas.jnc_counts_w_flags_&sysparm.;
   by fusion_id subject;
run;

proc means data=mysas.jnc_counts_w_flags_&sysparm. noprint;
   by fusion_id treatment;
   var apn log_apn region_depth;
   output out=jnc_counts_means mean=;
run;

*split on treatment, rename variables and merge;

data jnc_con_mean jnc_treat_mean;
    set jnc_counts_means;
    if treatment=0 then output jnc_con_mean;
    else if treatment=1 then output jnc_treat_mean;
    drop _TYPE_ _FREQ_ treatment;
run;

data jnc_con_mean_2;
   set jnc_con_mean;
   rename apn=mean_apn_con;
   rename log_apn=mean_log_apn_con;
   rename region_depth=mean_depth_con;
run;

data jnc_treat_mean_2;
   set jnc_treat_mean;
   rename apn=mean_apn_treat;
   rename log_apn=mean_log_apn_treat;
   rename region_depth=mean_depth_treat;
run;

data jnc_means_merge;
   merge jnc_con_mean_2 jnc_treat_mean_2;
   by fusion_id;
run;

/* making permenant */
/* adding in Up/down indicator for junctions. Doing this here instead of later for ease */
/* Can drop indicators later where junctions aren't diff expressed */

data mysas.jnc_means_merge_&sysparm.;
set jnc_means_merge;
   length flag_up_down $1.;
   if mean_apn_treat ge mean_apn_con then flag_up_down='U';
   else flag_up_down='D';
run;



