

proc freq data=fear.cis_data_estiamtes;
tables flag_ai_combined;
run;
*thse are all significant for ai;

proc freq data=fear.cis_data_estiamtes noprint;
tables fusion_id*line/out=count_fid_line;
run;

proc freq data=count_fid_line noprint;
tables fusion_id/out=count_line;
run;

proc sort data=fear.tests_anno;
by fusion_id;

proc sort data=count_line;
by fusion_id;

data check_counts;
merge fear.tests_anno count_line;
by fusion_id;
run;

proc freq data=check_counts;
tables count*flag_sig;
run;

proc freq data=check_counts;
where count>10 and flag_sig=1;
tables chrom*effect;
run;


proc freq data=check_counts;
where count>10 and effect="cis_i*trans_i";
tables chrom*flag_sig;
run;

data check_ago;
set check_counts;
where symbol_cat ? "AGO";
run;

data sig;
set check_counts;
where flag_all_ai=1;
rename count=num_lines;
run;


* no difference in X chromosome;
*genes on the x are not regulated dramatically differntly than those on the autosomes;

