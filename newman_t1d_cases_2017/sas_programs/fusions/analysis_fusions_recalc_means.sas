/* Recalculating means for analyzed data */

/* import libraries */
libname con '/home/jrbnewman/concannon/sas_data';

data cell_types;
  set con.design_file;

/* mhalanobis distance and low coverage*/
if name =  '2009-PC-0221' then delete; *sample 75 cd8;
if name =  '2009-PC-0144' then delete; *sample 48 cd4;
if name =  '2009-PC-0236' then delete; *sample 80;
if name =  '2009-PC-0237' then delete; *sample 80;
if name =  '2009-PC-0235' then delete; *sample 80;
 /*mahalnobis distance samples*/
  if name = '2009-PC-0101' then delete; *sample 34cd8;
if name = '2009-PC-0104' then delete; *sample35 cd8;
if name =  '2009-PC-0153' then delete; *sample 51 cd4;
if name =  '2009-PC-0200' then delete; *sample 67 cd8;
if name =  '2009-PC-0212' then delete; *sample 72 cd8;
if name = '2009-PC-0215' then delete; *sample 73  cd8;
/*funky heatmap samples*/

if name =  '2009-PC-0083' then delete; 
if name =   '2009-PC-0114' then delete;
if name =   '2009-PC-0224' then delete;
if name =   '2009-PC-0228' then delete;

  keep Name cell_type;
run;


proc sort data=con.fusion_q3_norm_data_all;
   by Name;
run;

proc sort data=cell_types;
by Name;
run;

data fusion_exp_w_cell oops1 oops2;
   merge cell_types (in=in1) con.fusion_q3_norm_data_all (in=in2);
   by Name;
   if in1 and in2 then output fusion_exp_w_cell;
   else if in1 then output oops1;
   else output oops2;
run;


proc sort data=fusion_exp_w_cell;
  by cell_type fusion_id;
run;

proc means data=fusion_exp_w_cell noprint;
   var log_q3_q3_apn;
   by cell_type fusion_id;
   output out=fusion_mean_counts_by_fusion mean=mean;
run;

/* Split means on cell_type */

data means_cd19 means_cd8 means_cd4 ;
   set fusion_mean_counts_by_fusion;
   if cell_type='CD19' then output means_cd19;
   if cell_type='CD4' then output means_cd4;
   if cell_type='CD8' then output means_cd8;
   drop _TYPE_ _FREQ_;
   run;

data means_cd19_2;
   set means_cd19;
   rename mean=mean_logq3q3apn_cd19;
   drop cell_type;
run;

data means_cd4_2;
   set means_cd4;
   rename mean=mean_logq3q3apn_cd4;
   drop cell_type;
run;

data means_cd8_2;
   set means_cd8;
   rename mean=mean_logq3q3apn_cd8;
   drop cell_type;
run;

data means_merge;
   merge means_cd19_2 means_cd4_2 means_cd8_2;
   by fusion_id;
run;


/* make permenant */

data con.fusion_means_by_cell_v2;
   set means_merge;
run;

