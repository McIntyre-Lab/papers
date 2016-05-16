libname cegs 'Z:/alison_g/cegs_ase_explore/sas_data';
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";

proc sort data=cegs.cis_calls;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
run;

data cis_calls_w_gene;
merge cegs.cis_calls (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run; *131798 obs;

*check how many genes have > 1 fusion id associated with them;
proc sort data=cis_calls_w_gene nodupkey out=by_fusion;
by fusion_id;
run; *3053 fusion_ids;

proc freq data=by_fusion;
tables symbol / out=chk_fusion_in_symbol;
run; *1593 genes;

data more_1_fusion;
set chk_fusion_in_symbol;
if count>1;
run; *755 have > 1 fusion ;

/****** OR DOES THIS WORK??? *********/
*this looks like what we want;
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

data all_calls_sbs;
merge all_v_sig_fusions2 all_m_sig_fusions2;
by line fusion_id;
run; *65948 obs;

proc sort data=all_calls_sbs ;
by line;
run;

proc freq data=all_calls_sbs noprint;
by line;
where flag_AI_M = 1 and flag_AI_V = 0;
tables symbol / out=count_M_enrich;
run;
proc freq data=all_calls_sbs noprint;
by line;
where flag_AI_V=1 and flag_AI_M = 0;
tables symbol / out=count_V_enrich;
run;
proc freq data=all_calls_sbs noprint;
by line;
where flag_AI_V = 1 and flag_AI_M = 1;
tables symbol / out=count_both_enrich;
run;
proc freq data=all_calls_sbs noprint;
by line;
where flag_AI_V=0 and flag_AI_M=0;
tables symbol / out=count_no_enrich;
run;

proc sort data=count_M_enrich;
by symbol;
proc sort data=count_V_enrich;
by symbol;
proc sort data=count_both_enrich;
by symbol;
run;

data gene_enrich_by_mating;
merge count_M_enrich (in=in1) count_V_enrich (in=in2) count_both_enrich (in=in3);
by symbol;
if in1 and in2 and in3 then gene_compare_enrich = "M_V_both";
else if in1 and in2 then gene_compare_enrich = "M_V";
else if in1 and in3 then gene_compare_enrich = "M_both";
else if in2 and in3 then gene_compare_enrich = "V_both";
else if in3 then gene_compare_enrich= "all_both";
else if in1 then gene_compare_enrich = "M_only";
else gene_compare_enrich= "V_only";
drop count percent;
run;

proc freq data=gene_enrich_by_mating;
tables gene_compare_enrich / out=chk_enrichments;
run; *by line and gene;
/* M_V			343
   M_V_both		6903
   M_both		293
   M_only		128
   V_both		398
   V_only		133
   all_both		293 */

proc freq data=gene_enrich_by_mating;
tables symbol*gene_compare_enrich / out=chk_enrich_gene;
run; *1605 obs;

proc sort data=chk_enrich_gene nodupkey out=unique_genes;
by symbol;
run; *0 obs deleted-- all genes have the same category of AI in all the lines;
* dataset unique_genes will have the gene and the category ;
*dataset gene_enrich_by_mating has each gene-line combination with a category;
*chk_enrich_gene is the same as unique_genes-- the count variable is the number of lines in the category;
*if the gene was in >1 category it would have had 2 lines so we checked that from the proc sort;
proc sort data=unique_genes;
by symbol;
proc sort data=all_calls_sbs;
by symbol;
run;

data unique_genes2;
merge unique_genes all_calls_sbs;
by symbol;
run;

data M_V_genes;
set unique_genes;
if gene_compare_enrich = "M_V";
run; *135 obs;

data M_V_both_genes;
set unique_genes;
if gene_compare_enrich = "M_V_both";
run; *913 obs;

data M_both_genes;
set unique_genes;
if gene_compare_enrich = "M_both";
run; *115 obs;

data M_only_genes;
set unique_genes;
if gene_compare_enrich = "M_only";
run; *93 obs;

data V_both_genes;
set unique_genes;
if gene_compare_enrich = "V_both";
run; *137 genes;

data V_only_genes;
set unique_genes;
if gene_compare_enrich = "V_only";
run; *99 obs;

data all_both_genes;
set unique_genes;
if gene_compare_enrich = "all_both";
run; *113 obs;

/* Gene numbers
   M_V			135
   M_V_both		913
   M_both		115
   M_only		93
   V_both		137
   V_only		99
   all_both		113 */

%macro exp (dat);
proc export data=&dat._genes
outfile = "Z:/cegs_ase_paper/output/gene_categories/&dat._genes.csv"
dbms=csv replace;
putnames=yes;
run;
%mend;
%exp (M_V);
%exp (M_V_both);		
%exp (M_both);	
%exp (M_only);	
%exp (V_both);	
%exp (V_only);	
%exp (all_both);

/* STOP */


/* below isn't quite what I want I think but I'm keeping it just in case */
data cis_mated;
set cis_calls_w_gene;
run;

/* must do by line first and then collapse onto gene */
proc sort data=cis_mated;
by line;
run;

proc freq data=cis_mated noprint;
by line;
tables mating_status*symbol/out=count_line_biased_genes;
run;

proc freq data=count_line_biased_genes noprint;
by line;
tables symbol /out=splice_count_line;
run;

data splice_diff_line;
set splice_count_line;
where count >1;
keep line symbol;
run;

*1379 with mixed bias;
proc sort data=count_line_biased_genes;
by line symbol;
proc sort data=splice_diff_line;
by line symbol;
data bias_count bias_diff oops;
merge count_line_biased_genes (in=in1) splice_diff_line (in=in2);
by line symbol;
if in1 and in2 then output bias_diff;
else if in1 then output bias_count;
else output oops;
run; *bias count=1075 obs, bias diff=82496 obs--- this is by LINE;

data sex_bias_line_1;
set bias_diff;
gene_bias_line="mixed";
keep line symbol gene_bias_line;

proc sort data=sex_bias_line_1 nodupkey out=sex_bias_line_1a;
by line symbol gene_bias_line;
run;

data sex_bias_line_2;
set bias_count;
gene_bias_line="bias_count";
keep line symbol gene_bias_line;
run;

data gene_bias_line_mated;
set sex_bias_line_2 sex_bias_line_1a ;
run;

proc sort data=gene_bias_line_mated;
by line symbol;
run;

proc freq data=gene_bias_line_mated;
tables gene_bias_line;
run; /* mated 
		bias count 1075
		mixed	   41248 */
*mixed means both mated and virgin share the genes, bias_count means the gene is particular to either mated or virgin;


proc sort data=gene_bias_line_mated;
by symbol;
proc sort data=cis_mated;
by symbol;
run;

data flag_gene_bias;
merge gene_bias_line_mated cis_mated;
by symbol;
if gene_bias_line = "bias_count" and flag_AI_combined=1 and mating_status="M" then flag_mated_only = 1;
else flag_mated_only=0;
if gene_bias_line= "bias_count" and flag_AI_combined=1 and mating_status="V" then flag_virgin_only=1;
else flag_virgin_only=0;
if flag_AI_combined =1 and gene_bias_line = "mixed" then flag_both =1;
else flag_both=0;
if flag_AI_combined = 0 then flag_neither = 1;
else flag_neither=0;
run;

*count unique mated only genes;
proc freq data=flag_gene_bias noprint;
where flag_mated_only=1 and flag_virgin_only=0;
tables symbol / out=count_mated_only;
run;

*count unique virgin only genes;
proc freq data=flag_gene_bias noprint;
where flag_virgin_only=1 and flag_mated_only=0;
tables symbol / out=count_virgin_only;
run;

proc freq data=flag_gene_bias noprint;
/* ??????????????????????????????????????????????????????? */


/* Now run process again collapsing on gene level by each flag */
data cis_mated_only;
set flag_gene_bias;
run;

proc freq data=cis_mated_only noprint;
tables mating_status*symbol/out=count_line_biased_mated_only;
run;

proc freq data=count_line_biased_mated_only noprint;
tables symbol /out=splice_count_mated_only;
run;

data splice_diff_mated_only;
set splice_count_mated_only;
where count >1;
keep symbol;
run;

*1379 with mixed bias;
proc sort data=count_line_biased_mated_only;
by symbol;
proc sort data=splice_diff_mated_only;
by symbol;
data bias_count_m_only bias_diff_m_only oops2;
merge count_line_biased_mated_only (in=in1) splice_diff_mated_only (in=in2);
by symbol;
if in1 and in2 then output bias_diff_m_only;
else if in1 then output bias_count_m_only;
else output oops2;
run; *bias count=207 obs, bias diff=2758 obs;

data sex_bias_line_m_only;
set bias_diff_m_only;
gene_bias_m_only="mixed";
keep symbol gene_bias_m_only;

proc sort data=sex_bias_line_m_only nodupkey out=sex_bias_m_only;
by symbol gene_bias_m_only;
run;

data sex_bias_line_m_only;
set bias_count_m_only;
gene_bias_m_only="bias_count";
keep symbol gene_bias_m_only;
run;

data gene_bias_line_m_only;
set sex_bias_line_m_only sex_bias_line_m_only ;
run;

proc sort data=gene_bias_line_m_only;
by symbol;
run;

proc freq data=gene_bias_line_m_only;
tables gene_bias_line2;
run; /* mated 
		bias count 187
		mixed	   1621 */
*mixed means both mated and virgin share the genes, bias_count means the gene is particular to either mated or virgin;

proc sort data=gene_bias_line_m_only;
by symbol;
proc sort data=flag_gene_bias;
by symbol;
run;

data flag_gene_bias2;
merge gene_bias_line_m_only flag_gene_bias;
by symbol;
if gene_bias_line2 = "bias_count" and flag_AI_combined=1 and mating_status="M" then flag_mated_only_gene = 1;
else flag_mated_only_gene=0;
if gene_bias_line2= "bias_count" and flag_AI_combined=1 and mating_status="V" then flag_virgin_only_gene=1;
else flag_virgin_only_gene=0;
if flag_AI_combined =1 and gene_bias_line2 = "mixed" then flag_both_gene =1;
else flag_both_gene=0;
if flag_AI_combined = 0 then flag_neither_gene = 1;
else flag_neither_gene=0;
run;

proc freq data=flag_gene_bias2;
where flag_mated_only=1 and flag_both=0 and flag_neither=0;
tables symbol / out=count_only_mated_enrich;
run;
proc freq data=flag_gene_bias2;
where flag_virgin_only=1 and flag_both=0 and flag_neither=0;


proc freq data=flag_gene_bias2;
tables flag_virgin_only_gene;
run;

proc freq data=flag_gene_bias2;
tables flag_virgin_only;
run;






proc freq data=flag_gene_bias;
tables flag_mated_only flag_virgin_only flag_both;
run;
/* flag_mated_only 0=25446, 1=189
flag_virgin_only 0=25512, 1=123
flag_both 0=312, 1=25323 */

data chk_mated;
set flag_gene_bias;
if flag_mated_only=1;
if q5_mean_theta > 0.5 then flag_line = 1;
else flag_line=0;
run;

proc freq data=chk_mated;
tables symbol*flag_line / out=chk_mated_gene_AI;
run; *145 obs or genes in mated with AI;

proc sort data=chk_mated_gene_AI nodupkey out=discordant_AI_mated_genes;
by symbol ;
run; *127 genes have AI in the same direction for each line represented;
*18 genes have discordant AI among the lines represetned;

proc freq data=cis_calls;
tables line*flag_AI_combined / out=chk_lines_w_AI;
run;

data percent_AI_r101;
set chk_lines_w_AI;
if line="r101" ;
run;
