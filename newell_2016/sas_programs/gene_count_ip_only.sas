libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";

proc freq data=ribo.gene_results_no_multi_with_exons;
  tables flag_ipfemale_on* flag_ipmale_on;
  run;
proc freq data=ribo.results_by_fusion_new_model;
  tables flag_ipfemale_on* flag_ipmale_on;
  run;

data results;
  set ribo.results_by_fusion_new_model;  where genes_per_fusion=1;

  run;

proc sort data=results nodupkey;
  by symbol_cat;
  run;
proc freq data=results;
  tables flag_ipfemale_on* flag_ipmale_on;
  run;




*trap differnt in two sexes;
proc freq data=all_results;
 where flag_fdr_p_contrast_1_20=1;
table sex_bias_ip*sex_bias;
run;

*306 exons shift  all are 3'most see the 3' program;
*ip only;
proc freq data=all_results noprint;
*significant sex differnce in ip;
where flag_fdr_p_contrast_1_20=1;
tables sex_bias_ip*symbol_cat/out=male_bias_ip_only_genes;
run;

data check1;
set male_bias_ip_only_genes;
where symbol_cat ? "fru";
run;
*fru is in this set;

proc freq data=male_bias_ip_only_genes noprint;
tables symbol_cat/out=count_sex_diff_ip_only;
run;

data splice_diff_ip_only;
set count_sex_diff_ip_only;
where count >1;
keep symbol_cat;
run;

*16 with mixed bias;

proc sort data=male_bias_ip_only_genes;
by symbol_cat;
proc sort data=splice_diff_ip_only;
by symbol_cat;

data bias_diff_ip_only bias_count_ip_only oops;
merge male_bias_ip_only_genes (in=in1) splice_diff_ip_only (in=in2);
by symbol_cat;
if in1 and in2 then output bias_diff_ip_only;
else if in1 then output bias_count_ip_only;
else output oops;
run;

data sex_bias_ip_only_1;
set bias_diff_ip_only;
gene_sexbias_ip_only="mixed";
keep symbol_cat gene_sexbias_ip_only;

proc sort data=sex_bias_ip_only_1 nodupkey out=sex_bias_ip_only_1a;
by symbol_cat gene_sexbias_ip_only;
run;

data sex_bias_ip_only_2;
set bias_count_ip_only;
gene_sexbias_ip_only=sex_bias_ip;
keep symbol_cat gene_sexbias_ip_only;
run;

data sex_bias_ip_only;
set sex_bias_ip_only_2 sex_bias_ip_only_1a;
run;


*getting 1741  total with uniform bias;
*16 bias diff;
* 0 oops!;
proc sort data=sex_bias_ip_only;
by symbol_cat;
proc freq data=sex_bias_ip_only;
tables gene_sexbias_ip_only;
run;

*1388 consitenly male biased in ip;
*353 consistently female biased in ip;
*16 mixed;

/* no multigene
female 331  
male 1210 
mixed 16 
*/
