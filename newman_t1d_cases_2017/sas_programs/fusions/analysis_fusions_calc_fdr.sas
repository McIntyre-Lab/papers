/* import libraries */
libname con '/home/jrbnewman/concannon/sas_data';

/* calculate FDR per comparison */



/* Prepare pairwise contrasts */

data cd4_8 cd4_19 cd8_19;
   set con.contrasts_fusions;
   if ProbF lt 0.05 then flag_p05=1;
   else flag_p05=0;
   if label='CD4-CD8' then output cd4_8;
   if label='CD4-CD19' then output cd4_19;
   if label='CD8-CD19' then output cd8_19;
   keep fusion_id ProbF flag_p05;
run;

data cd4_8_2;
   set cd4_8;
   rename ProbF=CD4_CD8_Pvalue;
   rename flag_p05=flag_p05_CD4_CD8;
run;

data cd4_19_2;
   set cd4_19;
   rename ProbF=CD4_CD19_Pvalue;
   rename flag_p05=flag_p05_CD4_CD19;
run;

data cd8_19_2;
   set cd8_19;
   rename ProbF=CD8_CD19_Pvalue;
   rename flag_p05=flag_p05_CD8_CD19;
run;


/* Calculate FDR */

proc multtest inpvalues(CD4_CD8_Pvalue)=cd4_8_2 fdr
 out=cd4_8_fdr noprint;
run;
quit;


data cd4_8_fdr2;
  set cd4_8_fdr;
  if fdr_p lt 0.05 then flag_cd4cd8_fdr05=1;
  else flag_cd4cd8_fdr05=0;
  rename fdr_p=CD4_CD8_FDR;
run; 

proc sort data=cd4_8_fdr2;
   by fusion_id;
run;

proc multtest inpvalues(CD4_CD19_Pvalue)=cd4_19_2 fdr
 out=cd4_19_fdr noprint;
run;
quit;

data cd4_19_fdr2;
  set cd4_19_fdr;
  if fdr_p lt 0.05 then flag_cd4cd19_fdr05=1;
  else flag_cd4cd19_fdr05=0;
  rename fdr_p=CD4_CD19_FDR;
run; 

proc sort data=cd4_19_fdr2;
   by fusion_id;
run;

proc multtest inpvalues(CD8_CD19_Pvalue)=cd8_19_2 fdr
 out=cd8_19_fdr noprint;
run;
quit;

data cd8_19_fdr2;
  set cd8_19_fdr;
  if fdr_p lt 0.05 then flag_cd8cd19_fdr05=1;
  else flag_cd8cd19_fdr05=0;
  rename fdr_p=CD8_CD19_FDR;
run; 

proc sort data=cd8_19_fdr2;
   by fusion_id;
run;

/* merge  */
data con.results_by_fusion_w_fdr oops1 oops2 oops3;
    merge cd4_8_fdr2 (in=in1) cd4_19_fdr2 (in=in2) cd8_19_fdr2 (in=in3);
    by fusion_id;
    if in1 and in2 and in3 then output con.results_by_fusion_w_fdr;
    else if in1 then output oops1;
    else if in2 then output oops2;
    else output oops3;
run;
