/* import libaries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* adding flags */

data counts_flag_on;
   set sugrue.counts_w_key;
   keep fusion_id treatment flag_fusion_on;
run;

proc sort data=counts_flag_on;
   by fusion_id treatment flag_fusion_on;
run;

proc freq data=counts_flag_on noprint;
   by fusion_id treatment;
   tables flag_fusion_on / out=counts_subjects_exp;
run;

*want a flag_on_control flag_on_treatment flag_on_all;

data counts_on_only;
   set counts_subjects_exp;
   if flag_fusion_on=1 then output;
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

data control_on;
    set control_on;
    if count>2 then flag_control_on=1;
    else flag_control_on=0;
    drop treatment count;
run;

data treat_on;
    set treat_on;
    if count>2 then flag_treat_on=1;
    else flag_treat_on=0;
    drop treatment count;
run;

proc sort data=control_on;
   by fusion_id;
run;

proc sort data=treat_on;
   by fusion_id;
run;

proc sort data=sugrue.counts_w_key;
   by fusion_id;
run;

data flags_on;
   merge control_on treat_on;
   by fusion_id;
   if flag_control_on=. then flag_control_on=0;
   if flag_treat_on=. then flag_treat_on=0;
   if flag_control_on=1 then do;
      if flag_treat_on=1 then flag_all_on=1;
      else flag_all_on=0;
      end;
   else flag_all_on=0;
run;

data counts_w_flags no_counts;
   merge sugrue.counts_w_key (in=in1) flags_on (in=in2);
   by fusion_id;
   if in1 and in2 then output counts_w_flags;
   else if in1 then do;
        flag_control_on=0;
        flag_treat_on=0;
        flag_all_on=0;
        output counts_w_flags;
        end;
   else output no_counts;
run;

/* Make permenant */ 

data sugrue.counts_w_flags;
   set counts_w_flags;
run;



