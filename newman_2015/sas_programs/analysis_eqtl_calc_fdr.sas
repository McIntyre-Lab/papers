/* Calc FDR for eQTLs */


libname eqtl '/mnt/data/eqtls/sas_data';

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


proc freq data=eqtl_results_w_fdr_flags;
    tables flag_eqtl_fdr05;
run;


/*
                                              Cumulative    Cumulative
  flag_eqtl_fdr05    Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------------
                0     1182995       99.12       1182995        99.12
                1       10558        0.88       1193553       100.00


*/

/* Make permenant */

data eqtl.eqtl_results_w_fdr;
   set eqtl_results_w_fdr_flags;
run;

