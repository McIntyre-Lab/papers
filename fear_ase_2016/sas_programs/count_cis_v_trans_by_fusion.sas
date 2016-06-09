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
out=standard_mating;
title2 'method=std';
var AI_effect;
by type;
run;

/***** mated cis *****/

*what's the range of cis-- flag low, med, neutral, high cis-- ;
*or flag if trans or cis is higher-- then compare proportion of trans vs cis higher in lines;
data cis_mating;
set stack_mating;
if type="cis_line";
rename AI_effect = cis_line;
run;

data trans_mating;
set stack_mating;
if type="trans_line";
rename AI_effect=trans_line;
run;

proc sort data=cis_mating;
by line fusion_id;
proc sort data=trans_mating;
by line fusion_id;
data mating_stdize;
merge cis_mating trans_mating;
by line fusion_id;
run;

data flag_cis_trans;
set mating_stdize;
if abs(cis_line) > abs(trans_line) then flag_cis_greater = 1;
else flag_cis_greater=0;
run;

proc freq data=flag_cis_trans;
tables fusion_id*flag_cis_greater / out=check_fusions;
run;

%macro prop (ID);
proc transpose data=check_fusions out=trans_fusion_&ID;
where fusion_id="&ID";
run;

data calc_prop_&ID;
set trans_fusion_&ID;
total_lines= COL1+COL2;
prop_cis = COL2/total_lines; *the proportion of lines that have cis>trans values;
fusion_id = "&ID";
prop_trans=1-prop_cis; *the proportion of lines that have trans>cis values;
keep prop_cis prop_trans fusion_id;
run;

data clean_calc_&ID;
set calc_prop_&ID (firstobs=2 obs=2);
run;
%mend;
%iterdataset(dataset=fusions, function=%nrstr(%prop(&fusion_id);));

data prop_cis_v_trans_fusions;
set clean_calc_: ;
run;

proc datasets noprint;
delete calc_prop_:;
delete trans_fusion_: ;
delete clean_calc_: ;
run; quit;

title 'Proportion lines per Fusion in AI for Mated Sigificant';
proc sgplot data=prop_cis_v_trans_fusions (obs=200);
where prop_cis > 0;
vbar fusion_id / response=prop_cis
categoryorder=respdesc;
run; quit; 

data flag_mating_cis;
length prop_cis_greater $ 20 ;
set prop_cis_v_trans_fusions;
if prop_cis = . then prop_cis_greater="no_diff";
else if prop_cis ge 0 and prop_cis le 0.2500 then prop_cis_greater = "low_le_025";
else if prop_cis ge 0.25100 and prop_cis le 0.500 then prop_cis_greater= "mid_le_05";
else if prop_cis ge 0.50100 and prop_cis le 0.7500 then prop_cis_greater = "high_le_075";
else if prop_cis ge 0.75100 then prop_cis_greater = "prevasive_ge_076";
run;

proc freq data= flag_mating_cis;
tables prop_cis_greater / out=cis_lines;
run;

/* %cis_greater #_fusions
	all 0 or 1	434
	0-25		2311
	26-50		281
	51-75		19
	>76			8  */

*this flag means that 0-25% of lines within 2311 fusions have cis>trans effects ;

data cegs.cnt_cis_trans_mated_by_fusion;
set flag_mating_cis;
run; *3035 obs;
