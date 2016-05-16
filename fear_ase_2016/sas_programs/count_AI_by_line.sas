filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/alison_g/cegs_ase_explore/sas_data';

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

*make cis_calls_sbs permanent;
data all_M_sig_fusions2;
set cis_calls_w_gene ;
if mating_status="M";
rename cis_line = cis_line_M;
rename trans_line = trans_line_M;
rename cis_tester = cis_tester_M;
rename trans_tester = trans_tester_M;
rename mean_apn= mean_apn_M;
rename flag_AI_combined = flag_AI_M;
rename q5_mean_theta = q5_mean_theta_M;
label cis_line = cis_line_M;
label trans_line = trans_line_M;
label cis_tester = cis_tester_M;
label trans_tester = trans_tester_M;
label mean_apn= mean_apn_M;
label flag_AI_combined = flag_AI_M;
label q5_mean_theta = q5_mean_theta_M;
drop mating_status mu;
run;

data all_V_sig_fusions2;
set cis_calls_w_gene ;
if mating_status="V";
rename cis_line = cis_line_V;
rename trans_line = trans_line_V;
rename cis_tester = cis_tester_V;
rename trans_tester = trans_tester_V;
rename mean_apn= mean_apn_V;
rename flag_AI_combined = flag_AI_V;
rename q5_mean_theta = q5_mean_theta_V;
label cis_line = cis_line_V;
label trans_line = trans_line_V;
label cis_tester = cis_tester_V;
label trans_tester = trans_tester_V;
label mean_apn= mean_apn_V;
label flag_AI_combined = flag_AI_V;
label q5_mean_theta = q5_mean_theta_V;
drop mating_status mu;
run;

proc sort data=all_v_sig_fusions2;
by line fusion_id;
proc sort data=all_m_sig_fusions2;
by line fusion_id;
run;

data cegs.all_calls_sbs;
merge all_v_sig_fusions2 all_m_sig_fusions2;
by line fusion_id;
run; *65948 obs;

data all_calls_sbs;
set cegs.all_calls_sbs;
run; *65948 obs;

proc sort data=all_calls_sbs ;
by line;
run;

data significant_AI;
set all_calls_sbs;
run; *18024 obs;

/*For mated*/
proc freq data=significant_AI;
tables line*flag_AI_M / out=chk_M_line;
run;

data lines;
set all_calls_sbs;
keep line;
run;

proc sort data=lines nodupkey;
by line;
run; *49 obs;

%macro prop (ID);
proc transpose data=chk_M_line out=trans_M_&ID;
where line="&ID";
run;

data calc_prop_&ID;
set trans_M_&ID;
total_regions= COL1+COL2;
prop_regions_AI_M = COL2/total_regions;
line = "&ID";
keep prop_regions_AI_M line;
run;

data clean_calc_&ID;
set calc_prop_&ID (firstobs=2 obs=2);
run;
%mend;
%iterdataset(dataset=lines, function=%nrstr(%prop(&line);));

data proportion_AI_lines;
set clean_calc_: ;
run;

proc datasets noprint;
delete calc_prop_:;
delete trans_M_: ;
delete clean_calc_: ;
run; quit;

/*For Virgin*/
proc freq data=significant_AI;
tables line*flag_AI_V / out=chk_V_line;
run;

%macro prop (ID);
proc transpose data=chk_V_line out=trans_V_&ID;
where line="&ID";
run;

data calc_prop_&ID;
set trans_V_&ID;
total_regions= COL1+COL2;
prop_regions_AI_V = COL2/total_regions;
line = "&ID";
keep prop_regions_AI_V line;
run;

data clean_calc_&ID;
set calc_prop_&ID (firstobs=2 obs=2);
run;
%mend;
%iterdataset(dataset=lines, function=%nrstr(%prop(&line);));

data proportion_AI_lines2;
set clean_calc_: ;
run;

proc datasets noprint;
delete calc_prop_:;
delete trans_M_: ;
delete clean_calc_: ;
run; quit;

proc sort data=proportion_AI_lines2;
by line;
proc sort data=proportion_AI_lines;
by line;
run;

data merge_prop;
merge proportion_AI_lines2 proportion_AI_lines;
by line;
run;

proc sort data=merge_prop;
by descending prop_regions_AI_V ;
run;

proc export data=merge_prop
outfile = "Z:/cegs_ase_paper/output/proportion_lines_w_AI.csv"
dbms=csv replace;
putnames=yes;
run;

proc transpose data=merge_prop out=trans_prop
name=proportion
label=proportion;
id line;
run;

data trans_prop2;
set trans_prop;
drop proportion;
run;

proc export data=trans_prop2
outfile = "Z:/cegs_ase_paper/output/proportion_lines_w_AI_trans.csv"
dbms=csv replace;
putnames=yes;
run;


title 'Proportion lines in AI for Mated Sigificant';
proc sgplot data=proportion_AI_lines;
vbar line / response=prop_regions_AI
categoryorder=respdesc;
run; quit;

data cegs.count_AI_mated_by_line;
length regions_with_AI $ 15;
set proportion_AI_lines;
if prop_regions_AI ge 0 and prop_regions_AI le 0.100 then regions_with_AI = "low_le_01";
else if prop_regions_AI ge 0.1001 and prop_regions_AI le 0.200 then regions_with_AI = "mid_le_02";
else if prop_regions_AI ge 0.2001 and prop_regions_AI le 0.300 then regions_with_AI ="high_le_03";
else if prop_regions_AI ge 0.3001 then regions_with_AI = "highest_ge_031";
run;

proc freq data=cegs.count_AI_mated_by_line;
tables regions_with_AI / out=count_AI_table;
run;
/* %regions_AI	#lines
	0-10		2
	11-20		28
	21-30		14
	>31			5  */

*by tester and line?????;
/* PLOT PROPORTION Q5 > 0.5 */
data significant_AI_M;
set all_calls_sbs;
if flag_AI_M = 1;
if q5_mean_theta_M > 0.5 then flag_to_line=1;
else flag_to_line=0;
run; *18024 obs;

proc freq data=significant_AI_M;
tables line*flag_to_line / out=chk_M_line2;
run;

data lines;
set all_calls_sbs;
keep line;
run;

proc sort data=lines nodupkey;
by line;
run; *49 obs;

%macro prop (ID);
proc transpose data=chk_M_line2 out=trans_M_&ID;
where line="&ID";
run;

data calc_prop_&ID;
set trans_M_&ID;
total_regions= COL1+COL2;
prop_regions_AI_g05 = COL2/total_regions;
line = "&ID";
keep prop_regions_AI_g05 line;
run;

data clean_calc_&ID;
set calc_prop_&ID (firstobs=2 obs=2);
run;
%mend;
%iterdataset(dataset=lines, function=%nrstr(%prop(&line);));

data proportion_AI_lines_g05;
set clean_calc_: ;
run;

proc datasets noprint;
delete calc_prop_:;
delete trans_M_: ;
delete clean_calc_: ;
run; quit;

title 'Proportion of Significant Mated q5 > 0.5';
proc sgplot data=proportion_AI_lines_g05;
vbar line / response=prop_regions_AI_g05
categoryorder=respdesc;
run; quit;

*compare the q>0.5 and proportion with AI;
proc sort data=proportion_ai_lines;
by line;
proc sort data=proportion_ai_lines_g05;
by line;
run;

data proportions_combined;
merge proportion_ai_lines proportion_ai_lines_g05;
by line;
run;

title 'Proportions of regions with AI vs Proportion of Sig Regions q>0.5 Mated';
proc gplot data=proportions_combined;
plot prop_regions_AI*prop_regions_AI_g05;
run;
quit;

proc corr data=proportions_combined;
var prop_regions_AI prop_regions_AI_g05;
run;
