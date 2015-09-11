/* import libaries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';


/* Calculate group means */

data data_for_means;
    set sugrue.counts_w_flags;
run;

proc sort data=data_for_means;
  by treatment fusion_id;
run;

proc means data=data_for_means noprint;
   var apn;
   by treatment fusion_id;
   output out=mean_counts_by_fusion sum=;
run;


data means_con means_treat;
   set mean_counts_by_fusion;
   if treatment=0 then output means_con;
   if treatment=1 then output means_treat;
   drop _TYPE_ _FREQ_;
   run;

data means_con2;
   set means_con;
   rename apn=mean_apn_control;
   drop treatment;
run;

data means_treat2;
   set means_treat;
   rename apn=mean_apn_treat;
   drop treatment;
run;

data means_merge;
   merge means_con2 means_treat2;
   by fusion_id;
run;

/* Add up-down indicator */

data means_up_down;
   set means_merge;
   length flag_up_down $1.;
   if mean_apn_control le mean_apn_treat then do;
        flag_up_down="U";
        if mean_apn_control=0 then fold_change=.;
        else fold_change=mean_apn_treat/mean_apn_control;
        end;
   else do;
        flag_up_down="D";
        if mean_apn_treat=0 then fold_change=.;
        else fold_change=-mean_apn_control/mean_apn_treat;
        end;
   if mean_apn_control<10 then flag_low_exp_con=1;
   else flag_low_exp_con=0;
   if mean_apn_treat<10 then flag_low_exp_treat=1;
   else flag_low_exp_treat=0;
   if mean_apn_control<10 and mean_apn_treat<10 then flag_low_exp_both=1;
   else flag_low_exp_both=0;
run;

/* add back into dataset and make permenant */

proc sort data=means_up_down;
   by fusion_id;
run;

proc sort data=data_for_means;
   by fusion_id;
run;

data sugrue.counts_w_means oops1 oops2;
   merge data_for_means (in=in1) means_up_down (in=in2);
   by fusion_id;
   if in1 and in2 then output sugrue.counts_w_means;
   else if in1 then output oops1;
   else output oops2;
run;

