ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";
libname con "!PATCON/sas_data";

/* Set up quant splicing test:
  model will be:
    by gene_id
    log_apn = cell_type fusion_id cell_type*fusion_id / htype=1;

(basic model, mostly for testing if model will work

I will take three gene for testing
*/

data genes4test;
   set event.hg19_genes_for_ds_test;
   keep gene_id fusion_id;
run;

data fus_counts;
   set con.apn0_q3_norm;
   if flag_low_coverage=1 then delete;
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

   keep name fusion_id log_q3_q3_apn_filter;

run;

data covars;
   set con.all_covariates;
   keep name cell_type subject_id sex p2 pool;
run;

proc sort data=genes4test;
    by fusion_id;
proc sort data=fus_counts;
   by fusion_id;
run;

data test_set;
   merge genes4test (in=in1) fus_counts (in=in2);
    by fusion_id;
   if in1 and in2;
run;

proc sort data=covars;
   by name;
proc sort data=test_set;
   by name;
run;

data test_set_w_covar;
  merge test_set (in=in1) covars (in=in2);
   by name;
  if in1 and in2;
run;

proc sort data=test_set_w_covar;
   by gene_id fusion_id cell_type;
run;


ods listing close;
proc glimmix data=test_set_w_covar;
   by gene_id;
   class cell_Type fusion_id sex pool subject_id;
   model log_q3_q3_apn_filter = sex cell_type fusion_id sex*cell_type cell_type*fusion_id pool p2 / htype=1;
    random pool;
    random resid / subject=subject_id;
   output out=resid resid=resid pred=pred student=stu;
lsmeans cell_type*fusion_id/pdiff;
ods output tests1=anova lsmeans=lsmeans;
run;
quit;

data flag_sig_int;
   set anova;
   where effect="cell_type*fusion_id";
   if ProbF<0.05 and ProbF >0 then flag_p05=1;
   else flag_p05=0;
run;

ods listing;
proc freq data=flag_sig_int;
   tables flag_p05;
run;

/*
                                       Cumulative    Cumulative
  flag_p05    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------
         0        1588       42.14          1588        42.14
         1        2180       57.86          3768       100.00

*/

proc sort data=flag_sig_int;
   by ProbF;
run;

/* Make these permenant so I can look at them later */

data event.t1d_quant_splice_anova;
   set anova;
run;

data event.t1d_quant_splice_resid;
   set resid;
run;

data event.t1d_quant_splice_lsmeans;
   set lsmeans;
run;

data event.t1d_quant_splice_flag_p05;
   set flag_sig_int;
run;



