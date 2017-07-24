ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";

/* For the T1D case-only data, for each quant splicing test (all 3 cell types, CD4/8, CD4/19, CD8/19, count genes with
   1. sig fusion*cell_type effects (P<0.05)
   2. failed residuals 

   Then, calc FDR on results and count significant genes */

data gene_w_multi;
   set hg19.hg19_aceview_fusions_si_info;
   where flag_multigene=1;
   keep gene_id;
run;


%macro flag_gene(cells);

data resid;
  set event.t1d_quant_splice_resid_&cells.;
run;

proc sort data=resid;
   by gene_id;
proc sort data=gene_w_multi nodup;
  by gene_id;
run;

data resid_flag_multigene;
   merge resid (in=in1) gene_w_multi (in=in2);
   by gene_id;
   if in2 then flag_multigene=1; else flag_multigene=0;
   if in1 then output;
run;

/* Flag and count residuals that failed normality -- fusion-level */
proc sort data=resid_flag_multigene;
   by fusion_id flag_multigene cell_type;
run;

proc univariate data=resid_flag_multigene normal noprint;
   by fusion_id flag_multigene ;
   var Resid;
   output out=normtest probn=pnorm;
run;

data flag_resids_fusion;
  set normtest;
  if pnorm = . then flag_fail_norm=.;
  else if pnorm le 0.05 then flag_fail_norm=1;
  else flag_fail_norm=0;
run;

proc freq data=flag_Resids_fusion ;
  tables flag_fail_norm;
run;

proc freq data=flag_Resids_fusion ;
  where flag_multigene=0;
  tables flag_fail_norm;
run;


/* Flag and count residuals that failed normality -- gene-level */
proc sort data=resid_flag_multigene;
   by gene_id flag_multigene fusion_id cell_type;
run;

proc univariate data=resid_flag_multigene normal noprint;
   by gene_id flag_multigene ;
   var Resid;
   output out=normtest probn=pnorm;
run;

data flag_resids_gene;
  set normtest;
  if pnorm = . then flag_fail_norm=.;
  else if pnorm le 0.05 then flag_fail_norm=1;
  else flag_fail_norm=0;
run;


proc freq data=flag_Resids_gene ;
  tables flag_fail_norm;
run;

proc freq data=flag_Resids_gene ;
  where flag_multigene=0;
  tables flag_fail_norm;
run;


/* Flag genes with significant fusion*cell interactions */

data flag_sig_int;
   set event.t1d_quant_splice_anova_&cells.;
   where effect="cell_type*fusion_id";
   if ProbF=. then flag_cell_by_fus_p05=.;
   else if ProbF<0.05 then flag_cell_by_fus_p05=1;
   else flag_cell_by_fus_p05=0;
run;

proc sort data=flag_sig_int;
  by gene_id;
proc sort data=gene_w_multi nodup;
  by gene_id;
run;

data flag_sig_int_w_mult;
  merge flag_sig_int (in=in1) gene_w_multi (in=in2);
  by gene_id;
  if in2 then flag_multigene=1; else flag_multigene=0;
  if in1 then output;
run;

proc freq data=flag_sig_int_w_mult;
   tables flag_cell_by_fus_p05;
run;

proc freq data=flag_sig_int_w_mult;
   where flag_multigene=0;
   tables flag_cell_by_fus_p05;
run;


/* Make these permenant */

data event.t1d_QS_flag_gene_p05_&cells.;
  set flag_sig_int_w_mult;
run;

data event.t1d_QS_flag_resid_fusion_&cells.;
  set flag_Resids_fusion;
run;

data event.t1d_QS_flag_resid_gene_&cells.;
  set flag_Resids_gene;
run;

%mend;

%flag_gene(all);
*%flag_gene(cd48);
*%flag_gene(cd419);
*%flag_gene(cd819);


/*

                                             Cumulative    Cumulative
  flag_fail_norm    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------------
               0         299        3.26           299         3.26
               1        8862       96.74          9161       100.00

                                             Cumulative    Cumulative
  flag_fail_norm    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------------
               0         158        4.19           158         4.19
               1        3610       95.81          3768       100.00


      flag_cell_                             Cumulative    Cumulative
      by_fus_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        3601       39.31          3601        39.31
               1        5560       60.69          9161       100.00

       flag_cell_                             Cumulative    Cumulative
       by_fus_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        1609       42.70          1609        42.70
                1        2159       57.30          3768       100.00
*/


