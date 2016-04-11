
proc freq data=all_results noprint;
where flag_fdr_p_contrast_5_20=1;
tables  sex_bias_ip*symbol_cat/out=count_sex_enrich_diff;
run;

data check;
set count_sex_enrich_diff;
where symbol_cat ? "fru";
run;
*fru is in this set;

proc freq data=count_sex_enrich_diff noprint;
tables symbol_cat/out=splice_count_sex_ip;
run;

data splice_diff_sex_ip;
set splice_count_sex_ip;
where count >1;
keep symbol_cat;
run;

*14 with mixed bias;

proc sort data=count_sex_enrich_diff;
by symbol_cat;
proc sort data=splice_diff_sex_ip;
by symbol_cat;

data bias_count_sex_ip bias_diff_sex_ip oops;
merge count_sex_enrich_diff (in=in1) splice_diff_sex_ip (in=in2);
by symbol_cat;
if in1 and in2 then output bias_diff_sex_ip;
else if in1 then output bias_count_sex_ip;
else output oops;
run;

data sex_bias_ip_1;
set bias_diff_sex_ip;
gene_sex_bias_ip="mixed";
keep symbol_cat gene_sex_bias_ip;

proc sort data=sex_bias_ip_1 nodupkey out=sex_bias_ip_1a;
by symbol_cat gene_sex_bias_ip;
run;

data sex_bias_ip_2;
set bias_count_sex_ip;
gene_sex_bias_ip=sex_bias_ip;
keep symbol_cat gene_sex_bias_ip;
run;

data sex_bias_ip_combined;
set sex_bias_ip_2 sex_bias_ip_1a ;
run;

proc sort data=sex_bias_ip_combined;
by symbol_cat;
run;

proc freq data=sex_bias_ip_combined;
tables gene_sex_bias_ip;
run;
*589 male ;
*163 female;
*14 mixed;

/*no multigene
female 151  
male 515  
mixed 14 */
