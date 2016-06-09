libname cegs "McLab/cegs_ase_paper/sas_data";
libname dmel "McLab/useful_dmel_data/flybase551/sasdata";

PROC IMPORT OUT= WORK.simulated_bias 
            DATAFILE= "McLab/cegs_ase_paper/manuscript/tables/sim_lines_bias_table.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC IMPORT OUT= WORK.qsim_bias 
            DATAFILE= "McLab/cegs_ase_paper/manuscript/tables/qsim_bias_table.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc sort data=simulated_bias;
by fusion_id;
proc sort data=qsim_bias;
by fusion_id;
run;

data full_table;
merge simulated_bias qsim_bias;
by fusion_id;
run;

*rename for clarification;
data full_table2;
set full_table;
rename fusion_id = Exonic_Region;
rename flag_exons_no_bias_simulated_lin = No_Bias_in_Simulation;
rename flag_exons_w_bias_simulated_line = Bias_in_Simulation;
rename flag_exons_bias_in_half_simulate = Bias_in_Half_of_Simulation;
rename flag_exons_bias_in_all_simulated = Bias_in_All_of_Simulation;
rename flag_exons_no_bias_qsim = No_Bias_in_qSim;
rename flag_exons_w_bias_qsim = Bias_in_qsim;
rename flag_exons_bias_in_half_qsim = Bias_in_Half_of_qsim;
rename flag_exons_bias_in_all_qsim = Bias_in_All_of_qsim;
run;

*convert all . to 0;
data nomiss(drop=i);                                                    
  set full_table2;                                                            
  array testmiss(*) _numeric_;                                            
  do i = 1 to dim(testmiss);                                              
    if testmiss(i)=. then testmiss(i)=0;                                    
  end;                                                                    
run;  

proc export data=nomiss
outfile = "McLab/cegs_ase_paper/manuscript/tables/Supplemental_Table1.csv"
dbms=csv replace;
putnames=yes;
run;

/*Check how many "only golden" fusions have had qSIM between (0.48-0.52)
how many of the remaining of these have the potential for bias from the genome ambiguity simulation */

data qsim_bayes;
set cegs.qsim_bayesian_results;
run;

data clean_stack;
set cegs.clean_ase_stack;
run;

proc freq data=clean_stack;
tables fusion_id*mating_status / out=chk_fusion;
run; *this is the "golden" data;

proc sort data=qsim_bayes;
by fusion_id line mating_status;
proc sort data=clean_stack;
by fusion_id line mating_status;
run;

data qsim_w_golden;
merge qsim_bayes clean_stack (in=in1);
by fusion_id line mating_status;
if in1;
run;

data flag_bias;
set qsim_w_golden;
if qsim_mean_theta > 0.48 and qsim_mean_theta < 0.52 then flag_no_bias =1;
else flag_no_bias=0;
run;

data flag_bias_m;
set flag_bias;
if mating_status="M";
run;

data flag_bias_v;
set flag_bias;
if mating_status="V";
run;

*for mated;
proc freq data=flag_bias_m;
tables fusion_id*flag_no_bias / out=chk_bias_m;
run;

data no_bias_m;
set flag_bias_m;
if flag_no_bias = 1;
run; *16258 obs;

proc freq data=no_bias_m;
tables fusion_id / out=chk_fusion_m;
run;

proc freq data=flag_bias_m;
tables fusion_id / out=chk_total_m;
run;

data chk_total_m2;
set chk_total_m;
rename COUNT = total_lines;
label COUNT = total_lines;
drop PERCENT;
run;

proc sort data=chk_total_m2;
by fusion_id;
proc sort data=chk_fusion_m;
by fusion_id;
run;

data compare;
merge chk_total_m2 chk_fusion_m;
by fusion_id;
if COUNT = . then COUNT = 0;
run;

data m_no_bias;
set compare;
percent_no_bias = COUNT/total_lines;
run;

data ten_per;
set m_no_bias;
if percent_no_bias > .10;
run; *3869 fusions had at least 10% of lines with no bias;

data twenty_per;
set m_no_bias;
if percent_no_bias > 0.20;
run; *2255 fusions had at least 20% of lines with no bias;

data fifty_per;
set m_no_bias;
if percent_no_bias > 0.50;
run; *114 fusions had at least 50% of lines with no bias;

*for virgin;
proc freq data=flag_bias_v;
tables fusion_id*flag_no_bias / out=chk_bias_v;
run;

data no_bias_v;
set flag_bias_v;
if flag_no_bias = 1;
run; *16258 obs;

proc freq data=no_bias_v;
tables fusion_id / out=chk_fusion_v;
run;

proc freq data=flag_bias_v;
tables fusion_id / out=chk_total_v;
run;

data chk_total_v2;
set chk_total_v;
rename COUNT = total_lines;
label COUNT = total_lines;
drop PERCENT;
run;

proc sort data=chk_total_v2;
by fusion_id;
proc sort data=chk_fusion_v;
by fusion_id;
run;

data compare;
merge chk_total_v2 chk_fusion_v;
by fusion_id;
if COUNT = . then COUNT = 0;
run;

data v_no_bias;
set compare;
percent_no_bias = COUNT/total_lines;
run;

data ten_per;
set v_no_bias;
if percent_no_bias > .10;
run; *3817 fusions had at least 10% of lines with no bias;

data twenty_per;
set v_no_bias;
if percent_no_bias > 0.20;
run; *2239 fusions had at least 20% of lines with no bias;

data fifty_per;
set v_no_bias;
if percent_no_bias > 0.50;
run; *94 fusions had at least 50% of lines with no bias;

data bias;
set flag_bias_v;
if flag_no_bias = 0;
run;

proc sort data=full_table;
by fusion_id;
proc sort data=bias;
by fusion_id;
run;

data check_bias;
merge full_table bias (in=in1);
by fusion_id;
if in1;
run;

/* combine cis and trans information with AI in exons, AI switch, and genes with splice AI for supplemental table */

PROC IMPORT OUT= WORK.exons_w_ai 
            DATAFILE= "McLab/cegs_ase_paper/manuscript/tables/exons_w_ai.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN; *mated and virgin side-by-side, 79,967 obs;
*NO-- use cegs.clean_ase_stack and merge on cis and trans counts;

data ai_switch;
set cegs.ai_switch;
run; *list of line, mating_status, symbol_cat. 365 obs;
*NO for supplemental;

data genes_splice_ai;
set cegs.genes_qith_spice_ai;
rename symbol_cat = Gene_Name;
rename num_lines_env = Number_of_Lines;
label num_lines_env = Number_of_Lines;
drop Percent;
run; *subset of ai_switch where it's just the symbol_cat or gene name ;

proc export data=genes_splice_ai
outfile ="McLab/cegs_ase_paper/output/splice_ai_gene_list.csv"
dbms=csv replace;
putnames=yes;
run;

/* Combine this list with the full results of AI with cis and trans
proc export data=genes_splice_ai
outfile = "McLab/cegs_ase_paper/manuscript/tables/Supplemental_Table4.csv"
dbms=csv replace;
putnames=yes;
run; */


data cis_trans;
set cegs.cis_est_v13;
keep fusion_id line mating_status c_i T_i_1a direction_cis direction_trans;
run; *mated and virgin stacked, 53,880 obs;

proc freq data=cis_trans;
tables direction_cis*direction_trans / out=chk_cis;
run;

proc export data=cis_trans
outfile = "McLab/cegs_ase_paper/output/cis_trans_estimates_v13.csv"
dbms=csv replace;
putnames=yes;
run;

data ai_results;
set cegs.clean_ase_stack;
keep line fusion_id mating_status q5_mean_theta flag_all_AI flag_AI_qsim flag_AI_combined sum_line sum_tester sum_both sum_total mean_apn;
run; *159934 obs;

proc sort data=cis_trans;
by fusion_id line mating_status;
proc sort data=ai_results;
by fusion_id line mating_status;
run;

data ai_w_cis_trans;
merge ai_results cis_trans;
by fusion_id line mating_status;
if c_i = . then flag_analyze_cis = 0;
else flag_analyze_cis = 1;
run; *159934 obs-- lots of missing cis and trans estimates... ;
*because cis and trans were only calculated if number of lines > 10;

data gene_info;
set dmel.Fb551_si_fusions_unique_flagged;
keep fusion_id chrom start end symbol_cat FBgn_cat ;
run;

proc sort data=ai_w_cis_trans;
by fusion_id;
proc sort data=gene_info;
by fusion_id;
run;

data w_genes;
merge ai_w_cis_trans (in=in1) gene_info;
by fusion_id;
if in1;
run;

data final_results;
retain line mating_status fusion_id q5_mean_theta flag_all_AI flag_AI_qsim flag_AI_combined sum_line sum_tester sum_both sum_total mean_apn flag_analyze_cis;
set w_genes;
rename line=genotype;
rename fusion_id = exonic_region;
rename q5_mean_theta= mean_AI;
rename flag_all_AI = AI_intersection_test;
rename flag_AI_qsim = AI_qsim_test;
rename flag_AI_combined = AI_intersection_qsim_both;
rename sum_line = counts_in_genotype;
rename sum_tester= counts_in_tester;
rename sum_both = counts_unassigned;
rename sum_total = counts_total;
rename c_i = cis_effects_genotype;
rename T_i_1a = trans_effects_genotype;
rename flag_analyze_cis = analyze_cis_effects;
rename chrom= chromosome;
rename start = start_site;
rename end= end_site;
rename symbol_cat = gene_name;
rename FBgn_cat = Flybase551_FBgn;
run;

*combine with genes_splice_ai or ai_switch;
data ai_switch2;
set ai_switch;
flag_AI_alternate_direction = 1;
rename symbol_cat= gene_name;
rename line = genotype;
run;

proc sort data=final_results;
by genotype mating_status gene_name;
proc sort data=ai_switch2;
by genotype mating_status gene_name;
run;

data final_results2;
merge final_results ai_switch2;
by genotype mating_status gene_name;
if flag_AI_alternate_direction ne 1 then flag_AI_alternate_direction =0;
run;

proc freq data=final_results2;
tables analyze_cis_effects*exonic_region / out=check_cis;
run; *rows 4516 to 5395= 879 exons;

proc export data=final_results2
outfile = "McLab/cegs_ase_paper/manuscript/tables/Supplemental_Table3.csv"
dbms=csv replace;
putnames=yes;
run;


/* Add mating vs virgin significance -- can't put all together because of multiple genes per fusion */
PROC IMPORT OUT= WORK.mating_genes 
            DATAFILE= "McLab/cegs_ase_paper/output/gene_categories/M_only_genes.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
PROC IMPORT OUT= WORK.virgin_genes 
            DATAFILE= "McLab/cegs_ase_paper/output/gene_categories/V_only_genes.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC IMPORT OUT= WORK.r2_cis_genes 
            DATAFILE= "McLab/cegs_ase_paper/output/gene_categories/top_r2_cis_genes.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
PROC IMPORT OUT= WORK.r2_trans_genes 
            DATAFILE= "McLab/cegs_ase_paper/output/gene_categories/top_r2_trans_genes.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
PROC IMPORT OUT= WORK.r2_int_genes 
            DATAFILE= "McLab/cegs_ase_paper/output/gene_categories/top_r2_interaction_genes.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data mating2;
set mating_genes;
significant_in_mated_only = 1;
rename symbol = gene_name;
drop var5 gene_compare_enrich count percent;
run;

data virgin2;
set virgin_genes;
significant_in_virgin_only = 1;
rename symbol = gene_name;
drop var5 gene_compare_enrich count percent;
run;

data cis2;
set r2_cis_genes;
if flag_M = 1 then high_cis_contribute_mated = 1;
else high_cis_contribute_mated = 0;
if flag_V = 1 then high_cis_contribute_virgin = 1;
else high_cis_contribute_virgin = 0;
rename symbol_cat = gene_name;
drop fusion_id flag_M flag_V flag_both;
run;

data trans2;
set r2_trans_genes;
if flag_M = 1 then high_trans_contribute_mated = 1;
else high_trans_contribute_mated = 0;
if flag_V = 1 then high_trans_contribute_virgin = 1;
else high_trans_contribute_virgin = 0;
rename symbol_cat = gene_name;
drop fusion_id flag_M flag_V flag_both var6;
run;

data int2;
set r2_int_genes;
if flag_M = 1 then high_interact_contribute_mated = 1;
else high_interact_contribute_mated = 0;
if flag_V = 1 then high_interact_contribute_virgin = 1;
else high_interact_contribute_virgin = 0;
rename symbol_cat = gene_name;
drop fusion_id flag_M flag_V flag_both var6;
run;

proc sort data=virgin2;
by gene_name ;
proc sort data=mating2;
by gene_name ;
proc sort data=cis2;
by gene_name;
proc sort data=trans2;
by gene_name;
proc sort data=int2;
by gene_name;
run;

*can't do all together with supp table 3-- so here this is a separate list;
data mating_virgin_genes;
length gene_name $ 25;
merge cis2 mating2 virgin2 trans2 int2;
by gene_name ;
if significant_in_virgin_only ne 1 then significant_in_virgin_only =0;
if significant_in_mated_only ne 1 then significant_in_mated_only =0;
if high_interact_contribute_virgin ne 1 then high_interact_contribute_virgin =0;
if high_interact_contribute_mated ne 1 then high_interact_contribute_mated=0;
if high_cis_contribute_virgin ne 1 then high_cis_contribute_virgin =0;
if high_cis_contribute_mated ne 1 then high_cis_contribute_mated =0;
if high_trans_contribute_virgin ne 1 then high_trans_contribute_virgin =0;
if high_trans_contribute_mated ne 1 then high_trans_contribute_mated =0;
run; *450 obs;

proc export data=mating_virgin_genes
outfile = "McLab/cegs_ase_paper/manuscript/tables/Supplemental_Table5.csv"
dbms=csv replace;
putnames=yes;
run;


data sig_exon;
set exons_w_ai;
if flag_AI_combined_m_or_v =1 ;
run; *11,803 obs;


/* count number of genes with the golden 5391 fusions */
data sbs;
set cegs.clean_ase_sbs;
run;

proc sort data=dmel.Fb551_si_fusions_unique_flagged;
by fusion_id;
proc sort data=sbs;
by fusion_id;
run;

data merge_gene;
merge sbs (in=in1) dmel.Fb551_si_fusions_unique_flagged(in=in2);
by fusion_id;
if in1;
run;

proc sort data=merge_gene out=no_dup nodupkey;
by fusion_id;
run;

proc freq data=no_dup;
tables genes_per_fusion / out=chk_num_genes;
run;
/* genes_per_fusion COUNT
	1	4963
	2	391
	3	31
	4	6 
So... 4963 + 2*391 + 3*31 + 4*6 = 5,862 genes represented in these fusions */

/*count number of exons significant for AI in each line */
data sig_ai;
set ai_results;
if flag_ai_combined = 1;
run;

proc freq data=sig_ai;
tables line / out=count_line;
run;

proc sort data=count_line;
by count;
run;
/*highest with AI = w47, 1811 exons
lowest with AI = r21, 173 exons*/

proc freq data=ai_results;
tables line / out=count_total_line;
run;

data count_total_line2;
set count_total_line;
rename count=total_lines;
label count=total_lines;
drop percent;
run;

proc sort data=count_total_line2;
by line;
proc sort data=count_line;
by line;
run;

data merge_line;
merge count_total_line2 count_line;
by line;
percent_in_AI = COUNT / total_lines;
run;

proc sort data=merge_line;
by percent_in_AI;
run;

proc export data=merge_line
outfile = "McLab/cegs_ase_paper/output/fusions_with_signficant_AI_summary.csv"
dbms=csv replace;
putnames=yes;
run;
