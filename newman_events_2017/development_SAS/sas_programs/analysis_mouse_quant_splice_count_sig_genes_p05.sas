ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";

/* For mouse NPC vs OLD quant splicing test, count genes with
   1. sig fusion*cell_type effects (P<0.05)
   2. failed residuals (these probably don't matter too much, since we have 2 reps for 2 cell types.
                        Can't *really* determine if residuals are normally distributed from 4 samples...

   Then, calc FDR on results and count significant genes */

data resid;
  set event.mouse_quant_splice_resid;
run;

/* Flag and count residuals that failed normality -- fusion-level */
proc sort data=resid;
   by fusion_id sample_id;
run;

proc univariate data=resid normal noprint;
   by fusion_id;
   var Resid;
   output out=normtest probn=pnorm;
run;

data flag_resids_fusion;
  set normtest;
  if pnorm = . then flag_fail_norm=.;
  else if pnorm le 0.05 then flag_fail_norm=1;
  else flag_fail_norm=0;
run;

proc freq data=flag_Resids_fusion noprint;
  tables flag_fail_norm / out=flag_fail_norm;
run;

proc print data=flag_fail_norm;
run;

/*
          flag_
          fail_
   Obs     norm    COUNT    PERCENT

    1       .          1      .
    2       0      68883    96.4640
    3       1       2525     3.5360

*/

/* Flag and count residuals that failed normality -- gene-level */
proc sort data=resid;
   by gene_id fusion_id sample_id;
run;

proc univariate data=resid normal noprint;
   by gene_id;
   var Resid;
   output out=normtest probn=pnorm;
run;

data flag_resids_gene;
  set normtest;
  if pnorm = . then flag_fail_norm=.;
  else if pnorm le 0.05 then flag_fail_norm=1;
  else flag_fail_norm=0;
run;

proc freq data=flag_Resids_gene noprint;
  tables flag_fail_norm / out=flag_fail_norm;
run;

proc print data=flag_fail_norm;
run;

/*
         flag_
         fail_
  Obs     norm    COUNT    PERCENT

   1       0       5251    83.7213
   2       1       1021    16.2787

*/


/* Flag genes with significant fusion*cell interactions */

data flag_sig_int;
   set event.mouse_quant_splice_anova;
   where effect="cell_type*fusion_id";
   if ProbF=. then flag_cell_by_fus_p05=.;
   else if ProbF<0.05 then flag_cell_by_fus_p05=1;
   else flag_cell_by_fus_p05=0;
run;

proc freq data=flag_sig_int;
   tables flag_cell_by_fus_p05;
run;


/*
    flag_cell_                             Cumulative    Cumulative
    by_fus_p05    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0        5703       90.93          5703        90.93
             1         569        9.07          6272       100.00
*/

/* Make these permenant */

data event.mm10_NPCvOLD_quantspl_flag_p05;
  set flag_sig_int;
run;

data event.mm10_NPCvOLD_quantspl_resid_gene;
  set flag_Resids_gene;
run;

data event.mm10_NPCvOLD_quantspl_resid_fus;
  set flag_Resids_fusion;
run;

