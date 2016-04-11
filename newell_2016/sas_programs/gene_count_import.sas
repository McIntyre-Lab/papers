
libname ribo 'C:\a1stuff\ribo';
libname ribo 'z:\arbeitman\arbeitman_ribotag\sas_data';

libname ribo '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/sas_data';

/* make a gene based file*/
proc contents data=ribo.results_by_fusion_new_model;
run;

data all_results;
set ribo.results_by_fusion_new_model;
length sex_bias_input $7 sex_bias_ip $7;
where flag_fusion_all_on0=1 and Genes_per_fusion=1;

if InputMale_mean_rpkm ge InputFeMale_mean_rpkm then sex_bias_input="male";
	else sex_bias_input="female";

if IPmale_mean_rpkm ge IPFemale_mean_rpkm then sex_bias_ip="male";
	else sex_bias_ip="female";
if sex_bias_input=sex_bias_ip then sex_bias=sex_bias_input;
	else sex_bias="shift";

diff_ip=IPmale_mean_rpkm - IPFemale_mean_rpkm;

If IPMale_mean_rpkm ge InputMale_mean_rpkm then male_bias=1;
	else male_bias=0;
If IPFeMale_mean_rpkm ge InputFeMale_mean_rpkm then female_bias=1;
	else female_bias=0;

if flag_fdr_p_contrast_2_20=1 and male_bias=1 then male_enrich=1;
	else male_enrich=0;
if flag_fdr_p_contrast_3_20=1 and female_bias=1 then female_enrich=1;
	else female_enrich=0;

if male_enrich=1 and female_enrich=1 and flag_fdr_p_contrast_1_20=1 then flag_trap_sex_diff=1;
	else flag_trap_sex_diff=0;

if male_enrich=1 and female_enrich=1 and flag_fdr_p_contrast_5_20=1 then flag_trap_sex_diff_enrich=1;
	else flag_trap_sex_diff_enrich=0;
run;



  proc contents data=alL_results;
  run;

proc freq data=all_results noprint;
tables symbol_cat/out=count_genes;
run;
*7889 genes total were analyzed;

