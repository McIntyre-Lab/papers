/*  We have a set of criteria for determining differentially expressed genes in the ribotag pull-down vs input data. This program will use these criteria to recreate the lists Michelle made. */

libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';
libname ribo '/home/fnew/mclab/arbeitman/arbeitman_ribotag/sas_data';
proc sort data=ribo.results_by_fusion;
  by fusion_id;
  run;

* (1) Only want fusions with genes_per_fusion=1, otherwise you cannot tell which overlapping gene is being expressed;
data one_gene_per_fusion;
  set ribo.results_by_fusion;
  if genes_per_fusion=1;
  run;

* (2) Background list, mean RPKM > 1 in at least one condition;
data background;
  set one_gene_per_fusion;
  if IPmale_mean_rpkm >1 or IPfemale_mean_rpkm >1 or InputMale_mean_rpkm >1 or InputFemale_mean_rpkm >1 then output background;
  run;
/* 46228 fusions have mean RPKM > 1 in at least one condition */

proc sort data=background nodupkey;
  by symbol_cat out=background_symbol;
  run; /* 11301 collapsed genes have mean RPKM > 1 in at least one condition */


/*  For the rest of the lists, flag_fdr_20 must be =1 for the comparison(s). 
The only common criteria for the following lists is gene_per_fusion=1. 
Do not use WORK.background for the rest of the lists. 
I need to add fold changes into the data for each comparison */


* (3) Enriched in TRAP not sex specific;
*	  >1 FC in Male_IP_Male_Input and >1 FC in FemaleIP_Female_Input ;
data fdr20_nonsex_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_20 eq 1 and flag_fdr_p_contrast_2_20 eq 1 then output fdr20_nonsex_trap;
  run;

data fdr20_ns_trap_fold;
  set fdr20_nonsex_trap;
  fold_male = ipmale_mean_rpkm / inputmale_mean_rpkm;
  fold_female = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  run;

data fdr_ns_trap_fc;
  set fdr20_ns_trap_fold;
  if fold_male >1 and fold_female >1 then output fdr_ns_trap_fc;
  run;
*2,100 fusions enriched in TRAP not sex specific ;

proc sort data=fdr_ns_trap_fc nodupkey;
  by fbgn_cat;
  run;
*1,642 genes enriched in TRAP not sex specific ;


* (4) Higher in female TRAP;
*	MaleIP_FemaleIP FC <1;

data fdr20_f_ip;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_1_20 eq 1 then output fdr20_f_ip;
  run;

data fdr20_f_ip_fold;
  set fdr20_f_ip;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_f_ip_fc;
  set fdr20_f_ip_fold;
  if fold_ip < 1 then output fdr_f_ip_fc;
  run;
*52 fusions higher in feamle TRAP;

proc sort data=fdr_f_ip_fc nodupkey;
  by fbgn_cat;
  run;
*40 collapsed genes higher in female TRAP;



* (5) Higher in Male TRAP;
*	MaleIP_FemaleIP FC >1 ;

data fdr_m_ip_fc;
  set fdr20_f_ip_fold;
  if fold_ip > 1 then output fdr_m_ip_fc;
  run;
*4,750 fusions higher in male TRAP;

proc sort data=fdr_m_ip_fc nodupkey;
  by fbgn_cat;
  run;
*3,125 collapsed genes higher in male TRAP;



* (6) Enriched in Female TRAP;
*	FemaleIP_FemaleInput FC >1 AND MaleIP_FemaleIP FC <1 ;

data fdr20_fem_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_20 eq 1 and flag_fdr_p_contrast_1_20 eq 1 then output fdr20_fem_trap;
  run;

data fdr20_fem_trap_fold;
  set fdr20_fem_trap;
  fold_fem = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_fem_trap_fc;
  set fdr20_fem_trap_fold;
  if fold_fem > 1 and fold_ip < 1 then output fdr_fem_trap_fc;
  run;
*19 fusions are enriched in female TRAP;

proc sort data=fdr_fem_trap_fc nodupkey;
  by fbgn_cat;
  run;
*19 collapsed genes are enriched in female TRAP;



* (7) Enriched in Male TRAP;
*	MaleIP_MaleInput FC >1 AND MaleIP_FemaleIP FC >1;

data fdr20_male_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_2_20 eq 1 and flag_fdr_p_contrast_1_20 eq 1 then output fdr20_male_trap;
  run;

data fdr20_male_trap_fold;
  set fdr20_male_trap;
  fold_male = ipmale_mean_rpkm / inputmale_mean_rpkm;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_male_trap_fc;
  set fdr20_male_trap_fold;
  if fold_male > 1 and fold_ip > 1 then output fdr_male_trap_fc;
  run;
*1,808 fusions are enriched in male TRAP;

proc sort data=fdr_male_trap_fc nodupkey;
  by fbgn_cat;
  run;
*1,528 genes enriched in male TRAP;


* (8) >2 FC Enriched in Female TRAP;
*	FemaleIP_FemaleInput FC >1 AND MaleIP_FemaleIP FC <0.5;

data fdr_fem_trap_fc_2;
  set fdr20_fem_trap_fold;
  if fold_fem >1 and fold_ip < 0.5 then output fdr_fem_trap_fc_2;
  run;
*16 fusions >2 FC enriched in Female TRAP;

proc sort data=fdr_fem_trap_fc_2 nodupkey;
  by fbgn_cat;
  run;
*16 collapsed genes >2 FC enriched in Female TRAP;



* (8) >2 FC Enriched in Male TRAP;
*	MaleIP_MaleInput FC >1 AND MaleIP_FemaleIP FC >2 ;

data fdr_male_trap_fc_2;
  set fdr20_male_trap_fold;
  if fold_male >1 and fold_ip > 2 then output fdr_male_trap_fc_2;
  run;
*281 fusions >2 FC enriched in male TRAP;

proc sort data=fdr_male_trap_fc_2 nodupkey;
  by fbgn_cat;
  run;
*270 genes are >2 FC enriched in male TRAP;




* (9) Higher in Female IP compared the Female Input;
*	FemaleIP_FemaleInput FC >1;

data fdr20_f_ip_in;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_20 eq 1  then output fdr20_f_ip_in;
  run;

data fdr20_f_ip_in_fold;
  set fdr20_f_ip_in;
  fold_f_ip = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  run;

data fdr20_f_ip_in_fc;
  set fdr20_f_ip_in_fold;
  if fold_f_ip > 1 then output fdr20_f_ip_in_fc;
  run;
* 3,069 fusions;

proc sort data=fdr20_f_ip_in_fc nodupkey;
  by fbgn_cat;
  run;
* 2,241 genes;


* (10) Higher in Male IP compared to Male Input;
*	MaleIP_MaleInput FC > 1;

data fdr20_m_ipin;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_2_20 = 1 then output fdr20_m_ipin;
  run;

data fdr20_m_ipin_fold;
  set fdr20_m_ipin;
  fold_m_ip = ipmale_mean_rpkm / inputmale_mean_rpkm;
  run;

data fdr_m_ipin_fc;
  set fdr20_m_ipin_fold;
  if fold_m_ip > 1 then output fdr_m_ipin_fc;
  run;
* 4,350 fusions are higher in Male IP than male input;

proc sort data=fdr_m_ipin_fc nodupkey;
  by fbgn_cat;
  run;
* 3,099 collapsed genes are higher in Male IP than male input;
