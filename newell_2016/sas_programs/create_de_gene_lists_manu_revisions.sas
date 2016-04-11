/* For the ribotag manuscript revisions,
>I am going to recreate Michelle's gene lists using fdr 0.1 and 0.05

>What are the fold changes?

*/


libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";

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
  if ipmale_mean_rpkm > 1 or IPfemale_mean_rpkm > 1 or inputmale_mean_rpkm > 1 or inputfemale_mean_rpkm > 1 
        then output background;
  run; *46,228 fusions have mean RPKM > 1 in at lest one condition;

proc sort data=background nodupkey;
  by symbol_cat;
  run; *11,301 collapsed genes have mean RPKM > 1 in at least one condition;

***************************************************************************************************
/*  For the rest of the lists, I am going to try FDR values of .1 and .05 for the comparison(s). 
The only common criteria for the following lists is gene_per_fusion=1. 
Do not use WORK.background for the rest of the lists. 
I need to add fold changes into the data for each comparison */

* (3) Enriched in TRAP not sex specific;
*       >1 FC in Male_IP_Male_Input and >1 FC in FemaleIP_Female_Input;

data fdr10_nonsex_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_10 eq 1 and flag_fdr_p_contrast_2_10 eq 1 then output fdr10_nonsex_trap;
  run;

data fdr10_ns_trap_fold;
  set fdr10_nonsex_trap;
  fold_male = ipmale_mean_rpkm / inputmale_mean_rpkm;
  fold_female = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  run;

data fdr_ns_trap_fc;
  set fdr10_ns_trap_fold;
  if fold_male > 1 and fold_female > 1 then output fdr_ns_trap_fc;
  run;
*2,073 fusions enriched in TRAP, not sex specific, at FDR 0.10;

proc sort data=fdr_ns_trap_fc nodupkey;
  by fbgn_cat;
  run;
*1,626 genes enriched in TRAP, not sex specific, at FDR 0.10;


* Now for FDR 0.05 ; 
data fdr5_nonsex_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_05 eq 1 and flag_fdr_p_contrast_2_05 eq 1 then output fdr5_nonsex_trap;
  run;

data fdr5_ns_trap_fold;
  set fdr5_nonsex_trap;
  fold_male = ipmale_mean_rpkm / inputmale_mean_rpkm;
  fold_female = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  run;

data fdr_ns_trap_fc;
  set fdr5_ns_trap_fold;
  if fold_male > 1 and fold_female > 1 then output fdr_ns_trap_fc;
  run;
*2,037 fusions enriched in TRAP, not sex specific, at FDR 0.10;

proc sort data=fdr_ns_trap_fc nodupkey;
  by fbgn_cat;
  run;
*1,607 genes enriched in TRAP, not sex specific, at FDR 0.10;

**********************************************************************
* (4) Higher in female TRAP;
*       MaleIP_FemaleIP  FC < 1;

data fdr10_f_ip;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_1_10 eq 1 then output fdr10_f_ip;
  run;

data fdr10_f_ip_fold;
  set fdr10_f_ip;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_f_ip_fc;
  set fdr10_f_ip_fold;
  if fold_ip < 1 then output fdr_f_ip_fc;
  run; * 19 fusions higher in female TRAP;

proc sort data=fdr_f_ip_fc nodupkey;
  by fbgn_cat;
  run; * 12 genes higher in female TRAP;


*Now for FDR 0.05;

data fdr5_f_ip;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_1_05 eq 1 then output fdr5_f_ip;
  run;

data fdr5_f_ip_fold;
  set fdr5_f_ip;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_f_ip_fc;
  set fdr5_f_ip_fold;
  if fold_ip < 1 then output fdr_f_ip_fc;
  run; * 14 fusions higher in female TRAP;

proc sort data=fdr_f_ip_fc nodupkey;
  by fbgn_cat;
  run; * 9 genes higher in female TRAP;

********************************************************************
* (5) Higher in Male TRAP;
*       MaleIP_FemaleIP FC > 1;

data fdr_m_ip_fc;
  set fdr10_f_ip_fold;
  if fold_ip > 1 then output fdr_m_ip_fc;
  run;
* 2,192 fusions higher in male TRAP;

proc sort data=fdr_m_ip_fc nodupkey;
  by fbgn_cat;
  run;
* 1,681 genes are higher in male TRAP;


* Now for FDR 0.05;
data fdr_m_ip_fc;
  set fdr5_f_ip_fold;
  if fold_ip > 1 then output fdr_m_ip_fc;
  run;
* 944 fusions higher in male TRAP;

proc sort data=fdr_m_ip_fc nodupkey;
  by fbgn_cat;
  run;
*791 genes higher in male TRAP;


****************************************************************************
* (6) Enriched in Female TRAP;
*   FemaleIP_FemaleInput FC > 1 AND MaleIP_FemaleIP FC < 1;

data fdr10_fem_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_10 eq 1 and flag_fdr_p_contrast_1_10 eq 1 then output fdr10_fem_trap;
  run;

data fdr10_fem_trap_fold;
  set fdr10_fem_trap;
  fold_fem = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_fem_trap_fc;
  set fdr10_fem_trap_fold;
  if fold_fem > 1 and fold_ip < 1 then output fdr_fem_trap_fc;
  run;
* 5 fusions higher in female TRAP than female input;

proc sort data= fdr_fem_Trap_fc nodupkey;
  by fbgn_cat;
  run;
* 5 genes higher in female TRAP than female input;

*Now for FDR 0.05;
data fdr5_fem_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_05 eq 1 and flag_fdr_p_contrast_1_05 eq 1 then output fdr5_fem_trap;
  run;

data fdr5_fem_trap_fold;
  set fdr5_fem_trap;
  fold_fem = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_fem_trap_fc;
  set fdr5_fem_trap_fold;
  if fold_fem > 1 and fold_ip < 1 then output fdr_fem_trap_fc;
  run;
* 4 fusions higher in female TRAP than female input;

proc sort data= fdr_fem_Trap_fc nodupkey;
  by fbgn_cat;
  run;
* 4 genes higher in female TRAP than female input at FDR 0.05;

******************************************************************************
*(7) Enriched in Male TRAP ;
*       MaleIP_MaleInput FC > 1 AND MaleIP_FemaleIP FC > 1;

data fdr10_male_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_2_10 eq 1 and flag_fdr_p_contrast_1_10 eq 1 then output fdr10_male_trap;
  run;

data fdr10_male_trap_fold;
  set fdr10_male_trap;
  fold_male = ipmale_mean_rpkm / inputmale_mean_rpkm;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_male_trap_fc;
  set fdr10_male_trap_fold;
  if fold_male > 1 and fold_ip > 1 then output fdr_male_trap_fc;
  run;
*917 fusions are enriched in male TRAP;

proc sort data=fdr_male_trap_fc nodupkey;
  by fbgn_cat;
  run;
*834 genes enriched in male TRAP;

* Now for FDR 0.05;

data fdr5_male_trap;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_2_05 eq 1 and flag_fdr_p_contrast_1_05 eq 1 then output fdr5_male_trap;
  run;

data fdr5_male_trap_fold;
  set fdr5_male_trap;
  fold_male = ipmale_mean_rpkm / inputmale_mean_rpkm;
  fold_ip = ipmale_mean_rpkm / ipfemale_mean_rpkm;
  run;

data fdr_male_trap_fc;
  set fdr5_male_trap_fold;
  if fold_male > 1 and fold_ip > 1 then output fdr_male_trap_fc;
  run;
*408 fusions are enriched in male TRAP;

proc sort data=fdr_male_trap_fc nodupkey;
  by fbgn_cat;
  run;
*388 genes enriched in male TRAP;

******************************************************************************************
* (8) FC >2 Enriched in Female TRAP;
*      FemaleIP_FemaleInput FC > 1 AND MaleIP_FemaleIP FC <0.5;

data fdr_fem_trap_fc_2;
  set fdr10_fem_trap_fold;
  if fold_fem > 1 and fold_ip < 0.5 then output fdr_fem_trap_fc_2;
  run;
* 5 fusions >2 FC enriched in Female TRAP;

proc sort data=fdr_fem_trap_fc_2 nodupkey;
  by fbgn_cat;
  run;
* 5 genes >2 FC enriched in Female TRAP;

*Now for FDR 0.05;

data fdr_fem_trap_fc_2;
  set fdr5_fem_trap_fold;
  if fold_fem > 1 and fold_ip < 0.5 then output fdr_fem_trap_fc_2;
  run;
* 4 fusions;

proc sort data=fdr_fem_trap_fc_2 nodupkey;
  by fbgn_cat;
  run;
* 4 genes;

*********************************************************************************************
* (9) >2 FC Enriched in Male TRAP;
*       MaleIP_MaleInput FC > 1 AND MaleIP_FemaleIP FC >2;

data fdr_male_trap_fc_2;
  set fdr10_male_trap_fold;
  if fold_male >1 and fold_ip > 2 then output fdr_male_trap_fc_2;
  run;
*192 fusions >2 FC enriched in male TRAP;

proc sort data=fdr_male_trap_fc_2 nodupkey;
  by fbgn_cat;
  run;
*187 genes are >2 FC enriched in male TRAP;


*Now for FDR 0.05;
data fdr_male_trap_fc_2;
  set fdr5_male_trap_fold;
  if fold_male >1 and fold_ip > 2 then output fdr_male_trap_fc_2;
  run;
*125 fusions >2 FC enriched in male TRAP;

proc sort data=fdr_male_trap_fc_2 nodupkey;
  by fbgn_cat;
  run;
*124 genes are >2 FC enriched in male TRAP;


**************************************************************************************************
* (10)  Higher in Female IP compared the Female Input;
*	FemaleIP_FemaleInput FC >1;

data fdr10_f_ip_in;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_10 eq 1  then output fdr10_f_ip_in;
  run;

data fdr10_f_ip_in_fold;
  set fdr10_f_ip_in;
  fold_f_ip = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  run;

data fdr10_f_ip_in_fc;
  set fdr10_f_ip_in_fold;
  if fold_f_ip > 1 then output fdr10_f_ip_in_fc;
  run;
* 3,042 fusions;

proc sort data=fdr10_f_ip_in_fc nodupkey;
  by fbgn_cat;
  run;
* 2,225 genes;


*Now for FDR 0.05;
data fdr5_f_ip_in;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_3_05 eq 1  then output fdr5_f_ip_in;
  run;

data fdr5_f_ip_in_fold;
  set fdr5_f_ip_in;
  fold_f_ip = ipfemale_mean_rpkm / inputfemale_mean_rpkm;
  run;

data fdr5_f_ip_in_fc;
  set fdr5_f_ip_in_fold;
  if fold_f_ip > 1 then output fdr5_f_ip_in_fc;
  run;
* 3,014 fusions;

proc sort data=fdr5_f_ip_in_fc nodupkey;
  by fbgn_cat;
  run;
* 2,211 genes;

****************************************************************************************************
* (11) Higher in Male IP compared to Male Input;
*	MaleIP_MaleInput FC > 1;

data fdr10_m_ipin;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_2_10 = 1 then output fdr10_m_ipin;
  run;

data fdr10_m_ipin_fold;
  set fdr10_m_ipin;
  fold_m_ip = ipmale_mean_rpkm / inputmale_mean_rpkm;
  run;

data fdr_m_ipin_fc;
  set fdr10_m_ipin_fold;
  if fold_m_ip > 1 then output fdr_m_ipin_fc;
  run;
* 4,323 fusions are higher in Male IP than male input;

proc sort data=fdr_m_ipin_fc nodupkey;
  by fbgn_cat;
  run;
* 3,087 collapsed genes are higher in Male IP than male input;


*Now for FDR 0.05;
data fdr5_m_ipin;
  set one_gene_per_fusion;
  if flag_fdr_p_contrast_2_05 = 1 then output fdr5_m_ipin;
  run;

data fdr5_m_ipin_fold;
  set fdr5_m_ipin;
  fold_m_ip = ipmale_mean_rpkm / inputmale_mean_rpkm;
  run;

data fdr_m_ipin_fc;
  set fdr5_m_ipin_fold;
  if fold_m_ip > 1 then output fdr_m_ipin_fc;
  run;
* 4,293 fusions are higher in Male IP than male input;

proc sort data=fdr_m_ipin_fc nodupkey;
  by fbgn_cat;
  run;
* 3,075 collapsed genes are higher in Male IP than male input;
