
*number of genes michelle counting;
proc freq data=all_results noprint;
where male_enrich=1 and female_enrich=1 and flag_fdr_p_contrast_1_20=1;
*contrast 1 is sex diff in ip;
tables symbol_cat/out=count_sex_trap;
run;

data enrich_in_both_sex_diff_in_ip;
set count_sex_trap;
enrich_in_both_sex_diff_in_ip=1;
keep symbol_cat enrich_in_both_sex_diff_in_ip;
run;



*71 genes;

/*no multi 
64 */
