filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/alison_g/cegs_ase_explore/sas_data';

*Is AI common? For a given fusion: lines with AI / total lines measured;

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
set all_calls_sbs;
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
run;

proc freq data=mating_fusions;
tables fusion_id*flag_AI_combined / out=check_fusions;
run;

%macro prop (ID);
proc transpose data=check_fusions out=trans_fusion_&ID;
where fusion_id="&ID";
run;

data calc_prop_&ID;
set trans_fusion_&ID;
total_lines= COL1+COL2;
prop_lines_AI = COL2/total_lines;
fusion_id = "&ID";
keep prop_lines_AI fusion_id;
run;

data clean_calc_&ID;
set calc_prop_&ID (firstobs=2 obs=2);
run;
%mend;
%iterdataset(dataset=fusions, function=%nrstr(%prop(&fusion_id);));

data proportion_AI_fusions;
set clean_calc_: ;
run;

proc datasets noprint;
delete calc_prop_:;
delete trans_fusion_: ;
delete clean_calc_: ;
run; quit;

title 'Proportion lines per Fusion in AI for Mated Sigificant';
proc sgplot data=proportion_AI_fusions (obs=200);
where prop_lines_AI > 0;
vbar fusion_id / response=prop_lines_AI
categoryorder=respdesc;
run; quit;

data no_diff_AI;
set proportion_AI_fusions;
if prop_lines_AI = .;
lines_with_AI = "no_diff";
run; *514 fusions have lines that are all AI=0 or AI=1;

data low_AI;
set proportion_AI_fusions;
if prop_lines_AI ge 0 and prop_lines_AI le 0.25;
lines_with_AI = "low_le_025";
run; *1835 fusions have between 0-25% of lines with AI;

data mid_AI;
set proportion_AI_fusions;
if prop_lines_AI ge 0.26 and prop_lines_AI le 0.5;
lines_with_AI= "mid_le_05";
run; *552 fusions have between 26-50% of lines with AI;

data high_AI;
set proportion_AI_fusions;
if prop_lines_AI ge 0.51 and prop_lines_AI le 0.75;
lines_with_AI = "high_le_075";
run; *111 fusions have between 51-75% of lines with AI;

data all_AI;
set proportion_AI_fusions;
if prop_lines_AI ge 0.76;
lines_with_AI = "prevasive_ge_076";
run; *23 fusions have AI with >76% of lines with AI;

/* %lines_w_AI #_fusions
	all 0 or 1	514
	0-25		1835
	26-50		552
	51-75		111
	>76			23  */

data cegs.count_AI_mated_by_fusion;
set all_AI High_AI Low_AI mid_AI no_diff_AI;
run; *3035 obs;
