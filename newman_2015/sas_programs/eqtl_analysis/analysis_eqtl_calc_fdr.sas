/* Check duplicate results */

libname eqtl '/mnt/data/eqtls/sas_data2';

data anova_results_for_fdr;
   set eqtl.anova_results;
   if effect='genotype';
run;


/* Quick check: are all feature*snp*cell_type entries unique? */

proc sort data=anova_results_for_fdr;
    by feature_id snp_id cell_type;
proc freq data=anova_results_for_fdr noprint;
    by feature_id snp_id;
    tables cell_type / out=eqtl_freq_check;
run;

data eqtl_dups;
   set eqtl_freq_check;
   if count gt 1;
run;


proc sort data=eqtl_dups;
    by feature_id snp_id cell_type;
run;

data dup_results;
   merge eqtl_dups (in=in1) anova_results_for_fdr;
   by feature_id snp_id cell_type;
   if in1;
run;


/* Okay, drop duplicated results ! */

proc sort data=anova_results_for_fdr nodup;
    by feature_id snp_id cell_type Fvalue ProbF;
run;


/* Make permenant */
data eqtl.anova_results_for_fdr;
   set anova_results_for_fdr;
run;

/* Calc FDR for eQTLs */
/* Parse anova results */

data eqtl_anova_results;
   set eqtl.anova_results_for_fdr;
   keep feature_id snp_id cell_type dendf ProbF;
run;

proc sort data=eqtl_anova_results;
  by feature_id snp_id cell_type probf;
run;


data eqtl_anova_results2;
   set eqtl_anova_results;
   by feature_id snp_id cell_type;
   if first.cell_type;
run;

proc sort data=eqtl_anova_results2;
   by descending dendf;
run;


/* Calculate FDR */

proc multtest inpvalues(ProbF)=eqtl_anova_results2 fdr
 out=eqtl_results_w_fdr noprint;
run;
quit;

data eqtl_results_w_fdr_flags;
   set eqtl_results_w_fdr;
   if ProbF=. then delete;
   else if ProbF lt 0.05 then flag_eqtl_p05=1;
   else flag_eqtl_p05=0;

   if fdr_p=. then delete;
   else do;
      if fdr_p lt 0.05 then flag_eqtl_fdr05=1;
      else flag_eqtl_fdr05=0;

      if fdr_p lt 0.20 then flag_eqtl_fdr20=1;
      else flag_eqtl_fdr20=0;
   end;
run;

ods listing; ods html close;
proc freq data=eqtl_results_w_fdr_flags;
    tables flag_eqtl_fdr05;
run;


/*

                                             Cumulative    Cumulative
 flag_eqtl_fdr05    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0     1193342       99.11       1193342        99.11
               1       10726        0.89       1204068       100.00

*/

/* Make permenant */

data eqtl.eqtl_results_w_fdr;
   set eqtl_results_w_fdr_flags;
run;

