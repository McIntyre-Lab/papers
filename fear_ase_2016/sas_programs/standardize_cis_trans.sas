filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/alison_g/cegs_ase_explore/sas_data';

filename mymacros "McLab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "McLab/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "McLab/useful_dmel_data/flybase551/sasdata";
libname cegs 'McLab/alison_g/cegs_ase_explore/sas_data';

*Is AI common? For a given fusion: lines with AI / total lines measured;

*proc standardize-- check in jmp-- sbs datafile-- normalize button;

proc sort data=cegs.cis_calls;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
run;

data cis_calls_w_gene;
merge cegs.cis_calls (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run;

data fusions;
set cis_calls_w_gene;
keep fusion_id;
run;

proc sort data=fusions nodupkey;
by fusion_id;
run; *3053 obs;

proc sort data=cis_calls_w_gene;
by line;
run;

/** MATED **/
data mating_fusions;
set cis_calls_w_gene;
if mating_status="M";
run; *65850 obs;

data cis_M;
set mating_fusions;
rename cis_line=AI_effect;
type="cis_line";
keep line mating_status fusion_id q5_mean_theta flag_AI_combined cis_line type;
run;

data trans_M;
set mating_fusions;
rename trans_line=AI_effect;
type="trans_line";
keep line mating_status fusion_id q5_mean_theta flag_AI_combined trans_line type;
run;

data stack_mating;
set trans_M cis_M;
run;

proc sort data=stack_mating;
by type;
run;

*visualize differences;
proc univariate data=stack_mating plot;
var AI_effect;
by type;
run;

*normalize using proc stdize and the mean and standard deviation as correction;
proc stdize data=stack_mating method=std pstat
out=cegs.standard_mating_stack;
title2 'method=std';
var AI_effect;
by type;
run;


/** VIRGIN **/
data virgin_fusions;
set cis_calls_w_gene;
if mating_status="V";
run; *65850 obs;

data cis_V;
set virgin_fusions;
rename cis_line=AI_effect;
type="cis_line";
keep line mating_status fusion_id q5_mean_theta flag_AI_combined cis_line type;
run;

data trans_V;
set virgin_fusions;
rename trans_line=AI_effect;
type="trans_line";
keep line mating_status fusion_id q5_mean_theta flag_AI_combined trans_line type;
run;

data stack_virgin;
set trans_V cis_V;
run;

proc sort data=stack_virgin;
by type;
run;

*visualize differences;
proc univariate data=stack_virgin plot;
var AI_effect;
by type;
run;

*normalize using proc stdize and the mean and standard deviation as correction;
proc stdize data=stack_virgin method=std pstat
out=cegs.standard_virgin_stack;
title2 'method=std';
var AI_effect;
by type;
run;
