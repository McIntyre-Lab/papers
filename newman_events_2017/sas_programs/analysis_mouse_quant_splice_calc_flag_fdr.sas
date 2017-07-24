ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";

/* For mouse NPC vs OLD quant splicing test, count genes with
   1. sig fusion*cell_type effects (P<0.05)
   2. failed residuals (these probably don't matter too much, since we have 2 reps for 2 cell types.
                        Can't *really* determine if residuals are normally distributed from 4 samples...

   Then, calc FDR on results and count significant genes */

data int_p05;
 set event.mm10_NPCvOLD_quantspl_flag_p05;
run;

proc multtest inpvalues(ProbF)=int_p05 fdr noprint
   out=int_p05_w_fdr;
run;

data flag_fdr05;
   set int_p05_w_fdr;
   if fdr_p=. then flag_cell_by_fus_fdr05=.;
   else if fdr_p<0.05 then flag_cell_by_fus_fdr05=1;
   else flag_cell_by_fus_fdr05=0;

   if fdr_p=. then flag_cell_by_fus_fdr10=.;
   else if fdr_p<0.10 then flag_cell_by_fus_fdr10=1;
   else flag_cell_by_fus_fdr10=0;


   if fdr_p=. then flag_cell_by_fus_fdr20=.;
   else if fdr_p<0.20 then flag_cell_by_fus_fdr20=1;
   else flag_cell_by_fus_fdr20=0;

run;

proc freq data=flag_fdr05;
  tables flag_cell_by_fus_fdr05 flag_cell_by_fus_fdr10 flag_cell_by_fus_fdr20;
run;

/*
   flag_cell_by_                             Cumulative    Cumulative
       fus_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        6119       97.56          6119        97.56
               1         153        2.44          6272       100.00


   flag_cell_by_                             Cumulative    Cumulative
       fus_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        6065       96.70          6065        96.70
               1         207        3.30          6272       100.00


   flag_cell_by_                             Cumulative    Cumulative
       fus_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        5974       95.25          5974        95.25
               1         298        4.75          6272       100.00
*/


/* Make FDR output permenant */
data event.mm10_NPCvOLD_quantspl_flag_fdr;
   set flag_fdr05;
run;


