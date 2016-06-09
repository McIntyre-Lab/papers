
libname cegs "S:\McIntyre_Lab\cegs_ase_paper\sas_data";

data cis_data_estimates;
set cegs.cis_data_estimates;
run;

proc sort data=cis_data_estimates;
by fusion_id mating_status;

ods noresults;
proc mixed data=cis_data_estimates  ;
by fusion_id mating_status ;
model bias=cis_i trans_i cis_i*trans_i/htype=1;
ods output tests1=tests1  ;
run;


data cegs.cis_trans_test;
set tests1;
if probf le 0 then flag_sig=.;
if probf >0 and probf<0.01 then flag_sig=1;
else flag_sig=0;
rename probf=raw_p;
run;

ods results;

proc freq data=cegs.cis_trans_test;
tables effect*flag_sig;
run;


proc multtest pdata=cegs.cis_trans_test fdr  out=fdr_cis;
run;

data cegs.fdr_flag_cis;
set fdr_cis;
if fdr_p le 0 then flag_fdr_sig=.;
if fdr_p >0 and fdr_p<0.05 then flag_fdr_sig=1;
else flag_fdr_sig=0;
run;

proc freq data=fdr_flag_cis;
tables effect*flag_fdr_sig;
run;
