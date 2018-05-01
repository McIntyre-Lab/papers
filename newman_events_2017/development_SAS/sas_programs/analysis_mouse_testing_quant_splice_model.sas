ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";

/* Set up quant splicing test:
  model will be:
    by gene_id
    log_apn = cell_type fusion_id cell_type*fusion_id / htype=1;

(basic model, mostly for testing if model will work

I will take three gene for testing
*/

data genes4test;
   set event.mm10_genes_for_ds_test;
   keep gene_id fusion_id;
run;

data fus_counts;
   length cell_type $3.;
   set event.mm10_refseq_fusion_counts;
   if sample_id = "NSC1" or sample_id = "NSC2" then cell_type="NSC";
   else cell_type="OLD";
   log_apn=log(apn+1);
   keep sample_id cell_type fusion_id apn log_apn;
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

proc sort data=test_set;
   by gene_id fusion_id cell_type;
run;

ods listing close;
proc glimmix data=test_set;
   by gene_id;
   class cell_Type fusion_id;
   model log_apn = cell_type|fusion_id / htype=1;
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
        0        5703       90.93          5703        90.93
        1         569        9.07          6272       100.00
*/

proc sort data=flag_sig_int;
   by ProbF;
run;

/* Make these permenant so I can look at them later */

data event.mouse_quant_splice_anova;
   set anova;
run;

data event.mouse_quant_splice_resid;
   set resid;
run;

data event.mouse_quant_splice_lsmeans;
   set lsmeans;
run;

data event.mouse_quant_splice_flag_p05;
   set flag_sig_int;
run;



