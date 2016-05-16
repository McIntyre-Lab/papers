libname cegs '!MCLAB/cegs_ase_paper/sas_data';

PROC IMPORT OUT= WORK.MBE_bayesian 
            DATAFILE= "/home/agerken/mclab/cegs_ase_paper/pipeline_output/emp_bayesian/MBE_PG_emp_bayesian_results.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data MBE_bayesian_M;
set MBE_bayesian;
if mating_status = 'M';
run;
*908679 obs;

data MBE_bayesian_V;
set MBE_bayesian;
if mating_status = 'V';
run;
*892330 obs;

PROC IMPORT OUT= WORK.counts_ase_bayes
            DATAFILE= "/home/agerken/mclab/cegs_ase_paper/pipeline_output/emp_bayesian/input/ase_dataset_for_bayesian.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data ase_bayesian_counts_M_all;
set counts_ase_bayes;
if mating_status = 'M';
if Line_total_4 = ' ' then num_reps = 3;
else num_reps = 6;
run;

data ase_bayesian_counts_M_3rep;
set ase_bayesian_counts_M_all;
if tester_total_4 = ' ' ;
Line_Average = (Line_total_1 + line_Total_2 + line_total_3 )/num_reps;
Tester_Average = (tester_total_1 + tester_Total_2 + tester_total_3)/num_reps;
Line_sum = Line_total_1 + line_Total_2 + line_total_3 ;
tester_sum = tester_total_1 + tester_Total_2 + tester_total_3;
Total_Count = Line_total_1 + line_Total_2 + line_total_3 + tester_total_1 + tester_Total_2 + tester_total_3;
Line_over_tester= Line_sum / Tester_sum;
Tester_over_line = Tester_sum / Line_sum;
run;

data ase_bayesian_counts_M_6rep;
set ase_bayesian_counts_M_all;
if tester_total_4 ne ' ' ;
Line_Average = (Line_total_1 + line_Total_2 + line_total_3 + line_total_4 + line_total_5 + line_total_6)/num_reps;
Tester_Average = (tester_total_1 + tester_Total_2 + tester_total_3 + tester_total_4 + tester_total_5 + tester_total_6)/num_reps;
Line_sum = Line_total_1 + line_Total_2 + line_total_3 + line_total_4 + line_total_5 + line_total_6;
Tester_sum = tester_total_1 + tester_Total_2 + tester_total_3 + tester_total_4 + tester_total_5 + tester_total_6;
Total_count = tester_total_1 + tester_Total_2 + tester_total_3 + tester_total_4 + tester_total_5 + tester_total_6 + Line_total_1 + line_Total_2 + line_total_3 + line_total_4 + line_total_5 + line_total_6;
Line_over_tester= Line_sum / Tester_sum;
Tester_over_line = Tester_sum / Line_sum;
run;

data ase_bayesian_counts_M;
set ase_bayesian_counts_M_6rep ase_bayesian_counts_M_3rep;
run;

data ase_bayesian_counts_V_all;
set counts_ase_bayes;
if mating_status = 'V';
if Line_total_4 = ' ' then num_reps = 3;
else num_reps = 6;
run;

data ase_bayesian_counts_V_3rep;
set ase_bayesian_counts_V_all;
if tester_total_4 = ' ' ;
Line_Average = (Line_total_1 + line_Total_2 + line_total_3 )/num_reps;
Tester_Average = (tester_total_1 + tester_Total_2 + tester_total_3)/num_reps;
Line_sum = Line_total_1 + line_Total_2 + line_total_3 ;
tester_sum = tester_total_1 + tester_Total_2 + tester_total_3;
Total_Count = Line_total_1 + line_Total_2 + line_total_3 + tester_total_1 + tester_Total_2 + tester_total_3;
Line_over_tester= Line_sum / Tester_sum;
Tester_over_line = Tester_sum / Line_sum;
run;

data ase_bayesian_counts_V_6rep;
set ase_bayesian_counts_V_all;
if tester_total_4 ne ' ' ;
Line_Average = (Line_total_1 + line_Total_2 + line_total_3 + line_total_4 + line_total_5 + line_total_6)/num_reps;
Tester_Average = (tester_total_1 + tester_Total_2 + tester_total_3 + tester_total_4 + tester_total_5 + tester_total_6)/num_reps;
Line_sum = Line_total_1 + line_Total_2 + line_total_3 + line_total_4 + line_total_5 + line_total_6;
Tester_sum = tester_total_1 + tester_Total_2 + tester_total_3 + tester_total_4 + tester_total_5 + tester_total_6;
Total_count = tester_total_1 + tester_Total_2 + tester_total_3 + tester_total_4 + tester_total_5 + tester_total_6 + Line_total_1 + line_Total_2 + line_total_3 + line_total_4 + line_total_5 + line_total_6;
Line_over_tester= Line_sum / Tester_sum;
Tester_over_line = Tester_sum / Line_sum;
run;

data ase_bayesian_counts_V;
set ase_bayesian_counts_V_6rep ase_bayesian_counts_V_3rep;
run;

proc sort data=MBE_bayesian_M;
by line fusion_id;
proc sort data=MBE_bayesian_V;
by line fusion_id;
proc sort data=ase_bayesian_counts_M;
by line fusion_id;
proc sort data=ase_bayesian_counts_V;
by line fusion_id;
run;

data MBE_bayes_M_w_counts;
merge MBE_bayesian_M (in=in1) ase_bayesian_counts_M;
by line fusion_id;
if in1;
run;
*908679 obs;

data MBE_bayes_V_w_counts;
merge MBE_bayesian_V (in=in1) ase_bayesian_counts_V;
by line fusion_id;
if in1;
run;
*892330 obs;

*for all together;
data MBE_bayes_all_w_counts;
set MBE_bayes_M_w_counts MBE_bayes_V_w_counts;
run;
*1801009 obs;

*make plots in R;
proc export data= MBE_bayes_all_w_counts 
outfile = '!MCLAB/cegs_ase_paper/r_data/emp_bayes_to_compare_20150218.csv'
dbms=csv replace;
putnames = yes;
run;

proc export data= MBE_bayes_M_w_counts 
outfile = '!MCLAB/cegs_ase_paper/r_data/emp_bayes_to_compare_M_20150218.csv'
dbms=csv replace;
putnames = yes;
run;

proc export data= MBE_bayes_V_w_counts 
outfile = '!MCLAB/cegs_ase_paper/r_data/emp_bayes_to_compare_V_20150218.csv'
dbms=csv replace;
putnames = yes;
run;



*for extended model;
PROC IMPORT OUT= WORK.counts_ase_bayes_extended
            DATAFILE= "/home/agerken/mclab/cegs_ase_paper/pipeline_output/extended_bayesian/input/data_for_bayes_BMC_PG_model.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

