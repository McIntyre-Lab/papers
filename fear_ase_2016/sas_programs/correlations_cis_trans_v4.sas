
libname cegs "S:\SHARE\\McIntyre_Lab\cegs_ase_paper\sas_data";
libname cegs '/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/sas_data';

*Correlation across fusions;


proc print data= cegs.cis_est_v13 (obs=10);
run;

data cis_est_mated;
set cegs.cis_est_v13;
where mating_status="M";
rename c_i=c_i_mated t_i_1a=t_i_1a_mated q5_mean_theta=q5_mean_theta_mated;
keep line fusion_id c_i t_i_1a q5_mean_theta;
run;

data cis_est_virgin;
set cegs.cis_est_v13;
where mating_status="V";
rename c_i=c_i_virgin t_i_1a=t_i_1a_virgin q5_mean_theta=q5_mean_theta_virgin;
keep line fusion_id c_i t_i_1a q5_mean_theta;
run;

data compare_est;
merge cis_est_mated (in=in1) cis_est_virgin(in=in2);
by fusion_id line;
if in1 and in2;
*assumption needed to compare environments;
run;

title compare raw estimates across environments;
proc gplot data= compare_est;
plot c_i_mated*c_i_virgin;
plot t_i_1a_mated*t_i_1a_virgin;
run;

proc sort data=compare_est;
by fusion_id;

proc corr data=compare_est out=all_by_fusion;
by fusion_id;
var c_i_mated c_i_virgin t_i_1a_mated t_i_1a_virgin;
run;



data corr_all_by_fusion;
set all_by_fusion;
where _type_="CORR";
run;

data cis_corr;
set corr_all_by_fusion;
where _NAME_="c_i_mated";
rename c_i_virgin=cis_corr;
keep fusion_id c_i_virgin;

data trans_corr;
set corr_all_by_fusion;
where _NAME_="t_i_1a_mated";
rename t_i_1a_virgin=trans_corr;
keep fusion_id t_i_1a_virgin;


data mated_corr;
set corr_all_by_fusion;
where _NAME_="c_i_mated";
rename t_i_1a_mated=ct_mated_corr;
keep fusion_id t_i_1a_mated;


data virgin_corr;
set corr_all_by_fusion;
where _NAME_="c_i_virgin";
rename t_i_1a_virgin=ct_virgin_corr;
keep fusion_id t_i_1a_virgin;


data corr_ct_by_e;
merge cis_corr trans_corr mated_corr virgin_corr;
by fusion_id;
if cis_corr>0 then direction_cis_corr="+";
else direction_cis_corr="-";
if trans_corr>0 then direction_trans_corr="+";
else direction_trans_corr="-";
abs_cis_corr=abs(cis_corr);
abs_trans_corr=abs(trans_corr);

if ct_mated_corr>0 then direction_ct_mated_corr="+";
else direction_ct_mated_corr="-";
if ct_virgin_corr>0 then direction_ct_virgin_corr="+";
else direction_ct_virgin_corr="-";

proc print data=corr_ct_by_e(obs=10);
run;

proc univariate data=corr_ct_by_e;
var ct_mated_corr ct_virgin_corr 
cis_corr trans_corr;
run;

proc gplot data=corr_ct_by_e;
plot cis_corr*trans_corr;
run;

proc freq data=corr_ct_by_e;
tables direction_cis_corr *direction_trans_corr;
tables direction_ct_mated_corr direction_ct_virgin_corr;
run;

data cegs.corr_ct_by_e;
set corr_ct_by_e;
run;

proc contents data=cegs.corr_ct_by_e;
run;

