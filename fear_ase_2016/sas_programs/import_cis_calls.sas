filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/cegs_ase_paper/sas_data';

data WORK.Cis_calls    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 'Z:/cegs_ase_paper/pipeline_output/cis_effects/cis_eq_full_output_w_ai_calls.csv' delimiter
  = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat line $4. ;
         informat mating_status $1. ;
         informat fusion_id $9. ;
         informat q5_mean_theta best32. ;
        informat flag_AI_combined best32. ;
         informat mean_apn best32. ;
         informat cis_tester best32. ;
        informat mu best32. ;
         informat cis_line best32. ;
         informat trans_tester best32. ;
        informat trans_line best32. ;
        format line $4. ;
        format mating_status $1. ;
        format fusion_id $9. ;
        format q5_mean_theta best12. ;
        format flag_AI_combined best12. ;
        format mean_apn best12. ;
         format cis_tester best12. ;
         format mu best12. ;
        format cis_line best12. ;
         format trans_tester best12. ;
        format trans_line best12. ;
      input
                  line $
                 mating_status $
                 fusion_id $
                  q5_mean_theta
                  flag_AI_combined
                 mean_apn
                  cis_tester
                  mu
                  cis_line
                  trans_tester
                  trans_line
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run; *131,700 records;

title 'cis_tester vs trans_tester';
proc gplot data=cis_calls;
plot cis_tester*trans_tester;
run;  quit; *these should line up because can't distinguish cis+trans in tester=0;

title 'cis_line vs trans_line';
proc gplot data=cis_calls;
plot cis_line*trans_line;
run; quit;

data sig_calls_M;
set cis_calls;
if flag_AI_combined=1 and mating_status="M";
run; *12653 obs;

title 'cis_line vs trans_line Mated Significant';
proc gplot data=sig_calls_M;
plot cis_line*trans_line;
run; quit;

data sig_calls_V;
set cis_calls;
if flag_AI_combined=1 and mating_status="V";
run; *12982 obs-- 
from JMF 491 regions with AI in mated and not virgin, 
469 regions with virgin and not mated;

title 'cis_line vs trans_line Virgin Significant';
proc gplot data=sig_calls_V;
plot cis_line*trans_line;
run; quit;

*check overlap in significant ??;

proc freq data=sig_calls_V;
tables fusion_id / out=chk_fusions_V;
run;

data many_fusions_V_sig;
set chk_fusions_V;
if COUNT > 10;
flag_V_sig = 1;
rename COUNT=V_COUNT;
run; *288 fusions have > 10 lines with significant AI ;

proc freq data=sig_calls_M;
tables fusion_id / out=chk_fusions_M;
run;

data many_fusions_M_sig;
set chk_fusions_M;
if COUNT > 10;
flag_M_sig = 1;
rename COUNT=M_Count;
run; *269 fusions have > 10 lines with significant AI;

proc sort data=many_fusions_M_sig;
by fusion_id;
proc sort data=many_fusions_V_sig;
by fusion_id;
run;

data merge_V_M;
merge many_fusions_M_sig many_fusions_V_sig;
if flag_M_sig ne 1 then flag_M_sig=0;
if flag_V_sig ne 1 then flag_V_sig=0;
run;

*iterdataset of the merge_v_m to see direction of cis and trans in these fusions;
*which is discordant between genotypes? Mated vs Virgin?
-- this is 2 questions first look at just M and just V, then look at M vs V;
*groups of lines?;
%macro make (fu);
data M_fusion_&fu;
set sig_calls_M;
if fusion_id = "&fu";
run;

%mend;
%iterdataset(dataset=many_fusions_M_sig, function=%nrstr(%make(&fusion_id);));

data all_M_sig_fusions;
set M_fusion_: ;
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

proc datasets noprint;
delete M_fusion_:;
run; quit;

%macro make (fu);
data V_fusion_&fu;
set sig_calls_V;
if fusion_id = "&fu";
run;

%mend;
%iterdataset(dataset=many_fusions_V_sig, function=%nrstr(%make(&fusion_id);));

data all_V_sig_fusions;
set V_fusion_: ;
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

proc datasets noprint;
delete V_fusion_:;
run; quit;

proc gplot data=all_m_sig_fusions;
plot cis_line*trans_line;
run; quit;

proc sort data=all_m_sig_fusions;
by line fusion_id;
proc sort data=all_v_sig_fusions;
by line fusion_id;
run;

data sig_fusions_m_v;
merge all_m_sig_fusions all_v_sig_fusions;
by line fusion_id;
if flag_AI_M + flag_AI_V > 1 then flag_AI_M_V = 1;
else flag_AI_M_V = 0;
run;

data sig_fusion_both_MV;
set sig_fusions_m_v;
if flag_AI_M_V = 1;
run;

proc sort data=all_M_sig_fusions;
by fusion_id;
proc sort data=all_V_sig_fusions;
by fusion_id;
run;

title 'mean apn vs cis line';
proc gplot data=cis_calls;
plot mean_apn*cis_line;
run; quit;

title 'mean apn vs cis line mated significant';
proc gplot data=all_m_sig_fusions;
plot mean_apn_M*cis_line_M;
run; quit;

*q5_mean_theta is direction of divergence to tester or line-- >0.5 to tester, <0.5 to line;
title 'q5_M vs q5_V significant AI';
proc gplot data=sig_fusion_both_mv;
plot q5_mean_theta_M*q5_mean_theta_V / vref=0.5 href=0.5;
run; quit; *there are several genotypes that are discordant where theta_M > theta_V and vice versa;

/* FIND FUSIONS THAT ARE DISCORDANT IN Q5 (BIAS) BETWEEN MATED AND VIRGIN */
data discordant;
set sig_fusion_both_mv;
if q5_mean_theta_M > 0.5 and q5_mean_theta_V < 0.5 then flag_MT_VL = 1;
else flag_MT_VL = 0;
if q5_mean_theta_M < 0.5 and q5_mean_theta_V > 0.5 then flag_ML_VT = 1;
else flag_ML_VT=0;
run;

data all_discrodant;
set discordant;
if flag_MT_VL =1 or flag_ML_VT=1;
run; *105 fusion-line discordant between mated and virgin;

data MT_VL;
set discordant;
if flag_MT_VL =1 ;
run; *61 fusion-line discordant between mated and virgin;

data ML_VT;
set discordant;
if flag_ML_VT=1;
run; *44 obs;

proc freq data=MT_VL;
tables fusion_id / out=chk_fusionid;
run; *very low number of lines in each fusion id;


proc sort data =all_discrodant;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
run;

data discordant_w_gene;
merge all_discrodant (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run;

/*** Export gene lists for significant in Mated only, significant in Virgin only, and overlap between Mated and Virgin **/

*how is cis distributed across the genome for significant for AI?? ;
proc sort data =all_m_sig_fusions;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
proc sort data=all_v_sig_fusions;
by fusion_id;
proc sort data=sig_fusion_both_MV;
by fusion_id;
run;

data m_sig_fusions_w_gene;
merge all_m_sig_fusions (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run; *4805 obs;

data v_sig_fusions_w_gene;
merge all_v_sig_fusions (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run;

data sig_both_MV_w_gene;
merge sig_fusion_both_MV (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run;

*export to plot-- this dataset is only the significant mated fusions that are unique to mated-- no overlap with virgin!;
proc export data= m_sig_fusions_w_gene
outfile = "Z:/cegs_ase_paper/output/m_sig_fusions_gene.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data= v_sig_fusions_w_gene
outfile = "Z:/cegs_ase_paper/output/v_sig_fusions_gene.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data= sig_both_MV_w_gene
outfile = "Z:/cegs_ase_paper/output/MV_overlap_sig_fusions_gene.csv"
dbms=csv replace;
putnames=yes;
run;

* now just by symbol and get an average for the gene of the theta;
proc sort data=sig_both_mv_w_gene;
by symbol;
proc means data=sig_both_mv_w_gene;
var q5_mean_theta_M ;
by symbol;
output out=m_overlap_mean mean=q5_mean_theta_M ;
run;

proc means data=sig_both_mv_w_gene;
var q5_mean_theta_V;
by symbol;
output out=v_overlap_mean mean=q5_mean_theta_V;
run;

data gene_list_mv_overlap;
merge v_overlap_mean m_overlap_mean;
by symbol;
run;

proc sort data=v_sig_fusions_w_gene;
by symbol;
proc means data=v_sig_fusions_w_gene;
var q5_mean_theta_v;
by symbol;
output out=gene_list_v_only mean=q5_mean_theta_V;
run;

proc sort data=m_sig_fusions_w_gene;
by symbol;
proc means data=m_sig_fusions_w_gene;
var q5_mean_theta_m;
by symbol;
output out=gene_list_m_only mean=q5_mean_theta_m;
run;

proc sort data=gene_list_mv_overlap;
by symbol;
proc sort data=dmel.genes2go_nodups;
by symbol;
proc sort data=gene_list_m_only;
by symbol;
proc sort data=gene_list_v_only;
by symbol;
run;

data m_only_w_go;
merge gene_list_m_only (in=in1) dmel.genes2go_nodups;
by symbol;
if in1;
run;
data v_only_w_go;
merge gene_list_v_only (in=in1) dmel.genes2go_nodups;
by symbol;
if in1;
run;
data mv_overlap_w_go;
merge gene_list_mv_overlap (in=in1) dmel.genes2go_nodups;
by symbol;
if in1;
run;

proc export data= mv_overlap_w_go
outfile = "Z:/cegs_ase_paper/output/gene_list_MV_overlap.csv"
dbms=csv replace;
putnames=yes;
run; *220 obs;

proc export data= m_only_w_go
outfile = "Z:/cegs_ase_paper/output/gene_list_M_only.csv"
dbms=csv replace;
putnames=yes;
run; *248 obs;

proc export data= v_only_w_go
outfile = "Z:/cegs_ase_paper/output/gene_list_V_only.csv"
dbms=csv replace;
putnames=yes;
run; *264 obs;




/* CHECK DISTRUTIONS AND COUNTS OF SIGNIFICANT FUSIONS */
proc freq data=m_sig_fusions_w_gene;
tables line / out=chk_line_m_sig;
run;

*this dataset is only fusions significant in M and not overlapping significant in V;
proc freq data=m_sig_fusions_w_gene;
where line="r101";
tables chrom / out=chk_chrom_r101;
run; /* chrom Count
		2L		18
		2R		31
		3L		15
		3R		9
		4		1
		X		15 */

proc freq data=m_sig_fusions_w_gene;
tables chrom / out=chk_chrom;
run; /* chrom Count
		2L		838
		2R		1713
		3L		883
		3R		584
		4		19
		X		768 */

/* CHECK DISTRIBUTION OF AI ACROSS THE CHROMOSOME-- DOES IT DIFFER IN CHROMOSOME LOCATION-- CHI SQUARED TESTS */
proc sort data=cis_calls;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
run;

data cis_calls_w_gene;
merge cis_calls (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run; *131798 obs;

proc sort data=cis_calls_w_gene;
by line;
run;
proc freq data=cis_calls_w_gene;
by line;
tables chrom*flag_AI_combined / out=chk_AI_dist_by_line;
run; 

proc freq data=cis_calls_w_gene;
by line;
tables chrom*flag_AI_combined / chisq;
run; /* to run in chi-sq test --should drop 4 and mitochondrion and 2R Het*/

data cis_calls_w_gene_clean;
set cis_calls_w_gene;
if chrom = "2RHet" or chrom= "4" or chrom = "dmel_mitochondrion_genome" then delete;
run; *went from 131798 to 88250 obs;

proc freq data=cis_calls_w_gene_clean;
by line;
tables chrom*flag_AI_combined / chisq;
ods output chiSq=chi_sq_output;
run; /* to run in chi-sq test without chrom 4, mitochondrion and 2RHet*/

data chi_sq_output_clean;
set chi_sq_output;
if statistic ne "Chi-Square" then delete;
run;

proc export data=chi_sq_output_clean
outfile="Z:/cegs_ase_paper/output/chi_sq_output_by_line.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=cis_calls_w_gene_clean
outfile="Z:/cegs_ase_paper/output/all_cis_calls_w_gene.csv"
dbms=csv replace;
putnames=yes;
run;

*check distribution of X vs not X;
data flag_X;
set cis_calls_w_gene_clean;
if chrom = "X" then flag_X=1;
else flag_X=0;
run;

proc sort data=flag_x;
by line;
run;

proc freq data=flag_x;
by line;
tables flag_X*flag_AI_combined / chisq;
ods output chiSq=chi_sq_output_X crossTabFreqs=Freq_Tables_X;
run; *typically there are a lot more AI events on not the X chromosome;

data chi_sq_out_X_clean;
set chi_sq_output_X;
if statistic ne "Chi-Square" then delete;
rename prob=chi_sq_pvalue;
keep line prob;
run;

data freq_tables_x_clean01;
set freq_tables_x;
where flag_X =0 and flag_AI_combined=1;
rename RowPercent = flagX_flagAI_01;
keep line RowPercent;
run;

data freq_tables_x_clean11;
set freq_tables_x;
where flag_X=1 and flag_AI_combined=1;
rename RowPercent = flagX_flagAI_11;
keep line RowPercent;
run;

proc sort data=freq_tables_x_clean01;
by line;
proc sort data=freq_tables_x_clean11;
by line;
proc sort data=chi_sq_out_x_clean;
by line;
run;

data freq_tables_x_clean_all;
merge freq_tables_x_clean11 freq_tables_x_clean01 chi_sq_out_x_clean;
by line;
diff_X_not_X = flagX_flagAI_11 - flagX_flagAI_01;
if diff_X_not_X > 0 then AI_higher_on_X =1;
else AI_higher_on_X = 0;
if chi_sq_pvalue < 0.05 then flag_sig_X_diff=1;
else flag_sig_x_diff=0;
run; *negative value means % with AI not on X is higher than % with AI on X;

proc univariate data=freq_tables_x_clean_all;
histogram diff_X_not_X;
run;

proc freq data=freq_tables_x_clean_all;
tables AI_higher_on_X*flag_sig_x_diff / out=check_AI_on_X;
run; /*AI_higher_on_X  flag_sig_x_diff Count (# lines)
		0				0				16
		1				1				10
		1				0				18
		1				1				5	*/
proc export data=freq_tables_x_clean_all
outfile="Z:/cegs_ase_paper/output/chi_sq_diff_AI_on_X.csv"
dbms=csv replace;
putnames=yes;
run;

*make file to check overlap location in M and V ASE;
data all_M_sig_fusions2;
set cis_calls_w_gene_clean ;
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
set cis_calls_w_gene_clean ;
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

data all_calls_sbs;
merge all_v_sig_fusions2 all_m_sig_fusions2;
by line fusion_id;
if flag_AI_V=1 and flag_AI_M = 1 then flag_AI_both=1;
else flag_AI_both=0;
run;

proc export data=all_calls_sbs
outfile="Z:/cegs_ase_paper/output/all_cis_calls_sbs.csv"
dbms=csv replace;
putnames=yes;
run;

