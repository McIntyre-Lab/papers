/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

/* Add in flags for counts etc. */

data counts_flag_on;
   set mysas.jnc_counts_w_key_&sysparm.;
   keep fusion_id treatment flag_junction_on;
run;

proc sort data=counts_flag_on;
   by fusion_id treatment flag_junction_on;
run;

proc freq data=counts_flag_on noprint;
   by fusion_id treatment;
   tables flag_junction_on / out=counts_subjects_exp;
run;

*want a flag_on_control flag_on_treatment flag_on_all;

data counts_on_only;
   set counts_subjects_exp;
   if flag_junction_on=1 then output;
   drop percent;
run;

proc sort data=counts_on_only;
  by fusion_id;
run;

data control_on treat_on;
   set counts_on_only;
   if treatment ne 1 then output control_on;
   if treatment ne 0 then output treat_on;
run;

data control_on_2;
    set control_on;
    if count>2 then flag_control_on=1;
    else flag_control_on=0;
    drop treatment count;
run;

data treat_on_2;
    set treat_on;
    if count>2 then flag_treat_on=1;
    else flag_treat_on=0;
    drop treatment count;
run;

proc sort data=control_on_2;
   by fusion_id;
run;

proc sort data=treat_on_2;
   by fusion_id;
run;

proc sort data=mysas.jnc_counts_w_key_&sysparm.;
   by fusion_id;
run;

data flags_on;
   merge control_on_2 treat_on_2;
   by fusion_id;
   if flag_control_on=. then flag_control_on=0;
   if flag_treat_on=. then flag_treat_on=0;
   if flag_control_on=1 then do;
      if flag_treat_on=1 then flag_all_on=1;
      else flag_all_on=0;
      end;
   else flag_all_on=0;
run;


data mysas.jnc_counts_w_flags_&sysparm. no_counts;
   merge mysas.jnc_counts_w_key_&sysparm. (in=in1) flags_on (in=in2);
   by fusion_id;
   if in1 and in2 then output mysas.jnc_counts_w_flags_&sysparm.;
   else if in1 then do;
        flag_control_on=0;
        flag_treat_on=0;
        flag_all_on=0;
        output mysas.jnc_counts_w_flags_&sysparm.;
        end;
   else output no_counts;
run;


/* flags for 10_counts */

data junc_10counts_flag;
   set mysas.jnc_counts_w_key_&sysparm.;
   if region_depth ge 10 then do;
        flag_10_cnts=1;
        flag_low_count=0;
        end;
   else do;
       flag_10_cnts=0;
       if region_depth gt 0 and region_depth lt 10 then flag_low_count=1;
       else flag_low_count=0;
       end;
   keep fusion_id treatment flag_10_cnts flag_low_count;
run;

proc freq data=junc_10counts_flag;
   tables flag_10_cnts flag_low_count;
run;



proc sort data=junc_10counts_flag;
   by fusion_id flag_10_cnts;
run;

proc freq data=junc_10counts_flag noprint;
   by fusion_id;
   tables flag_10_cnts / out=flag_10counts_freq;
run;

proc sort data=junc_10counts_flag;
   by fusion_id flag_low_count;
run;

proc freq data=junc_10counts_flag noprint;
   by fusion_id;
   tables flag_low_count / out=flag_lowcounts_freq;
run;


* pull junctions with 10counts or lowcounts;

data ten_counts_on_only;
   set flag_10counts_freq;
   if flag_10_cnts=1 then output;
   drop percent flag_10_cnts;
   rename count=num_samples_exp;
run;

data lowcounts_on_only;
   set flag_lowcounts_freq;
   if flag_low_count=1 then output;
   drop percent flag_low_count;
   rename count=num_samples_present;
run;

proc sort data=ten_counts_on_only;
  by fusion_id;
run;

proc sort data=lowcounts_on_only;
  by fusion_id;
run;

data junc_counts_merge;
    merge ten_counts_on_only lowcounts_on_only;
    by fusion_id;
    if num_samples_exp=. then num_samples_exp=0;
    if num_samples_present=. then num_samples_present=0;
run;

/* make permenant */

data mysas.junc_counts_merge_&sysparm.; *this is small. saving to share;
 set junc_counts_merge;
run;

