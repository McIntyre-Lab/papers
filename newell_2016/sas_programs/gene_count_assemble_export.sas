
libname sem 's:\McIntyre_Lab\cegs_sem_sd_paper\sas_data';

proc sort data=all_results;
by symbol_cat;

data gene_results;
merge all_results
gene_sex_bias_input
gene_enrich_by_sex 
sex_bias_ip_combined 
sex_bias_ip_only
enrich_in_both_sex_diff_in_ip ;
by symbol_cat;
run;

data gene_results2;
set gene_results;
if Gene_sex_bias_input="" then Gene_sex_bias_input="ns";
if Gene_compare_enrichment_bysex="" then Gene_compare_enrichment_bysex="ns";
if gene_sex_bias_ip="" then gene_sex_bias_ip="ns";
if gene_sexbias_ip_only="" then gene_sexbias_ip_only="ns";
if enrich_in_both_sex_diff_in_ip=. then enrich_in_both_sex_diff_in_ip=0;
run;

proc freq data=gene_results2;
tables Gene_sex_bias_input
Gene_compare_enrichment_bysex
gene_sex_bias_ip
gene_sexbias_ip_only
enrich_in_both_sex_diff_in_ip;
run;


McIntyre_Lab/cegs_sem_sd_paper/sas_data/cegsv_ag_yp2_flag_ds_fru_bic12.sas7bdat
flag_all_fru: 1 indicates at least one model was ds of fru
flag_best_fru: 1 indicates that best model was ds of fru

libname sem '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/cegs_sem_sd_paper/sas_data';

data ds_fru;
  set sem.cegsv_ag_yp2_flag_ds_fru_bic12;
  rename flag_all_fru= flag_downstream_fru symbol=symbol_cat;
  drop flag_best_fru primary_fbgn;
run;


proc sort data=ds_fru;
  by symbol_cat;
  run;



data gene_results_fru;
  merge gene_results2 (in=in1) ds_fru (in=in2);
  by symbol_cat;
  if in1;
  if flag_downstream_fru=. then flag_downstream_fru=0;
run;


proc freq data=gene_results_fru;
tables flag_downstream_fru*(Gene_sex_bias_input
Gene_compare_enrichment_bysex
gene_sex_bias_ip
gene_sexbias_ip_only
enrich_in_both_sex_diff_in_ip)/chisq;
run;

proc export data=gene_results_fru
	outfile="C:\a1stuff\ribo\gene_results_no_multi_with_exons.csv"
	DBMS=CSV REPLACE;
	run;

