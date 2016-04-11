
title male versus female enrichment in input;

proc freq data=all_results noprint;
*significant sex differnce in input;
where flag_fdr_p_contrast_4_20=1;
tables sex_bias_input*symbol_cat/out=count_male_bias_input_genes;
run;

proc freq data=count_male_bias_input_genes noprint;
tables symbol_cat/out=splice_count_input;
run;

data splice_diff_input;
set splice_count_input;
where count >1;
keep symbol_cat;
run;

*22 with mixed bias;
proc sort data=count_male_bias_input_genes;
by symbol_cat;
proc sort data=splice_diff_input;
by symbol_cat;
data bias_count bias_diff oops;
merge count_male_bias_input_genes (in=in1) splice_diff_input (in=in2);
by symbol_cat;
if in1 and in2 then output bias_diff;
else if in1 then output bias_count;
else output oops;
run;

data sex_bias_input_1;
set bias_diff;
gene_sex_bias_input="mixed";
keep symbol_cat gene_sex_bias_input;

proc sort data=sex_bias_input_1 nodupkey out=sex_bias_input_1a;
by symbol_cat gene_sex_bias_input;
run;

data sex_bias_input_2;
set bias_count;
gene_sex_bias_input=sex_bias_input;
keep symbol_cat gene_sex_bias_input;
run;

data gene_sex_bias_input;
set sex_bias_input_2 sex_bias_input_1a ;
run;

proc sort data=gene_sex_bias_input;
by symbol_cat;
run;

proc freq data=gene_sex_bias_input;
tables gene_sex_bias_input;
run;

*1685 male;
*1321 female;
*22 mixed;

/*no multigene
female 1125 43.55 1125 43.55 
male 1438 55.67 2563 99.23 
mixed 20 0.77 2583 100.00 */

