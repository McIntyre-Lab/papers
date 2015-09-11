/* import libraries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* Calculate FDR */
/* for only fusions on in both, and fusions on in either */


data fus_p_for_fdr;
    set sugrue.results_by_fusion;
    keep fusion_id flag_all_on ProbF;
run;


proc multtest inpvalues(ProbF)=fus_p_for_fdr fdr
 out=fusion_fdr_exp noprint;
  where flag_all_on=1;
run;
quit;


proc sort data=fusion_fdr_exp;
   by fusion_id;
run;

/* merge  */
data results_by_fusion_w_fdr oops1 oops2;
    merge sugrue.results_by_fusion (in=in1) fusion_fdr_exp (in=in2);
    by fusion_id;
    if in1 and in2 then output results_by_fusion_w_fdr;
    else if in1 then output results_by_fusion_w_fdr;
    else output oops2;
run;

/* add FDR flags and make permenant */
data sugrue.results_by_fusion_w_fdr2;
   set results_by_fusion_w_fdr;
   if fdr_p=. then flag_fdr_05=.;
   else if fdr_p<0.05 then flag_fdr_05=1;
   else flag_fdr_05=0;
run;



