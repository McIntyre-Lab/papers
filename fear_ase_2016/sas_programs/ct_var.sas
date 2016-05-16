 libname cegs "S:\SHARE\\McIntyre_Lab\cegs_ase_paper\sas_data";

 proc sort data=cegs.cis_est_v13;
 by fusion_id mating_status;
run;

 proc means data=cegs.cis_est_v13 noprint;
 by fusion_id mating_status;
 var c_i ;
 output out=cis mean=cis_mean var=cis_var;

proc means data=cegs.cis_est_v13 noprint;
 by fusion_id mating_status;
 var t_i_1a ;
 output out=trans mean=trans_mean var=trans_var;

 data compare_var;
 merge cis trans;
by fusion_id mating_status;

proc univariate data=compare_var normal plot;
var cis_var trans_var;
run;

proc print data=compare_var(obs=10);
run;

data cegs.compare_ct_var;
set compare_var;
run;

