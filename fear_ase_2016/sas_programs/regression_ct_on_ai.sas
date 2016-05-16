
libname cegs '/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/sas_data';


Data cis_int;
set cegs.cis_est_v13 ;
int= c_i*t_i_1a;

proc sort data=cis_int;
by fusion_id mating_status;

ods listing close;

*full model;

proc reg data=cis_int ;
by fusion_id  mating_status;
model q5_mean_theta=c_i t_i_1a int;
ods output ParameterEstimates=parms_full fitstatistics=fit_full;
run;

*no interaction;
proc reg data=cis_int ;
by fusion_id  mating_status;
model q5_mean_theta=c_i t_i_1a ;
ods output ParameterEstimates=parms_noint fitstatistics=fit_noint;
run;
quit;

*cis_only;

proc reg data=cis_int ;
by fusion_id  mating_status;
model q5_mean_theta=c_i  ;
ods output ParameterEstimates=parms_cis fitstatistics=fit_cis;
run;
quit;

data cegs.ai_reg_parms_full;
set parms_full;
if estimate ge 0 then sign_parm=”+”;
else sign_parm=”-”;

data cegs.ai_reg_parms_ct;
set parms_noint;
if estimate ge 0 then sign_parm=”+”;
else sign_parm=”-”;

data cegs.ai_reg_parms_cis;
set parms_cis;
if estimate ge 0 then sign_parm=”+”;
else sign_parm=”-”;

data cegs.ai_reg_fit_full2;
set cegs.ai_reg_fit_full;
where label2 ? ”Adj R-Sq”;
R2_full=cvalue2*1;
keep fusion_id mating_status R2_full;

data cegs.ai_reg_fit_noint2;
set cegs.ai_reg_fit_noint;
where label2 ? ”Adj R-Sq”;
R2_noint=cvalue2*1;
keep fusion_id mating_status R2_noint;
run;

data cegs.ai_reg_fit_cis;
set fit_cis;
where label2 ? ”Adj R-Sq”;
R2_cis=cvalue2*1;
keep fusion_id mating_status R2_cis;
run;

data cegs.r2_ct_models;
merge cegs.ai_reg_fit_full2 cegs.ai_reg_fit_noint2 cegs.ai_reg_fit_cis;
by fusion_id mating_status;
R2_diff_int=R2_full-R2_noint;
R2_diff_trans=R2_noint-R2_cis;
run;

ods listing;
proc univariate data =r2 normal plot  ;
var R2_full R2_noint R2_diff_int R2_diff_trans r2_cis;
run;



ods listing close;
proc reg data=cis_int ;
by fusion_id  mating_status;
model c_i=t_i_1a  ;
ods output ParameterEstimates=parms_cvt fitstatistics=fit_cvt;
run;
quit;

ods listing;

proc print data=cegs.parms_cvt (obs=10);
run;

data cegs.parms_cvt;
set parms_cvt;
where variable=”T_i_1a”;
run;

proc univariate data=cegs.parms_cvt;
var estimate;
run;
*at least 75% of estiamtes are negative; 

proc means data=r2 ;
by mating_status;
var  R2_full R2_noint R2_diff_int R2_diff_trans r2_cis;
run;

