ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";

/*  For the T1D case-only data, for each quant splicing test (all 3 cell types, CD4/8, CD4/19, CD8/19, calc FDR on results and count significant genes  */

%macro calc_fdr(cells);

data int_p05;
 set event.t1d_QS_flag_gene_p05_&cells.;
run;

proc multtest inpvalues(ProbF)=int_p05 fdr noprint
   out=int_p05_w_fdr;
   where flag_multigene=0;
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


/* Make FDR output permenant */

data event.t1d_QS_flag_gene_fdr_&cells.;
   set flag_fdr05;
run;

%mend;

%calc_fdr(all);
*%calc_fdr(cd48);
*%calc_fdr(cd419);
*%calc_fdr(cd819);

/*


      flag_cell_by_                             Cumulative    Cumulative
          fus_fdr05    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        1735       46.05          1735        46.05
                  1        2033       53.95          3768       100.00


      flag_cell_by_                             Cumulative    Cumulative
          fus_fdr10    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        1576       41.83          1576        41.83
                  1        2192       58.17          3768       100.00


      flag_cell_by_                             Cumulative    Cumulative
          fus_fdr20    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        1362       36.15          1362        36.15
                  1        2406       63.85          3768       100.00




*/

