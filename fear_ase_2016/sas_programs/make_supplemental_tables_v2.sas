libname cegs "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/sas_data";
libname dmel "/home/mcintyre/S/SHARE/McIntyre_Lab/useful_dmel_data/flybase551/sasdata";

PROC IMPORT OUT= WORK.simulated_bias 
            DATAFILE= "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/manuscript/tables/sim_lines_bias_table.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC IMPORT OUT= WORK.qsim_bias 
            DATAFILE= "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/manuscript/tables/qsim_bias_table.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc print data=dmel.Fb551_si_fusions_unique_flagged (obs=10);
run;


proc sort data=simulated_bias;
by fusion_id;
proc sort data=qsim_bias;
by fusion_id;
run;

data full_table;
merge simulated_bias qsim_bias;
by fusion_id;
run;

proc sort data=dmel.Fb551_si_fusions_unique_flagged;
by fusion_id;
run;

proc sort data=full_table;
by fusion_id;
run;

data simulation_results oops;
merge full_table (in=in1) dmel.Fb551_si_fusions_unique_flagged(in=in2);
by fusion_id;
if in1 and in2 then output simulation_results;
else output oops;
rename fusion_id=exonic_region;
run;


proc export data=simulation_results
outfile = "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/manuscript/resubmission/Supplemental_Table1_lmm2.csv"
dbms=csv ;
putnames=yes;
run;


/* Add mating vs virgin significance -- can't put all together because of multiple genes per fusion */
PROC IMPORT OUT= WORK.mating_genes 
            DATAFILE= "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/output/gene_categories/M_only_genes.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
PROC IMPORT OUT= WORK.virgin_genes 
            DATAFILE= "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/output/gene_categories/V_only_genes.csv" 
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

proc sort data=virgin2;
by gene_name ;
proc sort data=mating2;
by gene_name ;

*can't do all together with supp table 3-- so here this is a separate list;
data mating_virgin_genes;
length gene_name $ 25;
merge mating2 virgin2 ;
by gene_name ;
if significant_in_virgin_only ne 1 then significant_in_virgin_only =0;
if significant_in_mated_only ne 1 then significant_in_mated_only =0;
run; *191 obs;

proc export data=mating_virgin_genes
outfile = "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/manuscript/resubmission/Supplemental_Table4_lmm2.csv"
dbms=csv replace;
putnames=yes;
run;

proc sort data=cegs.r2_ct_models;
by fusion_id mating_status;
run;


data regression_results ;
merge  cegs.corr_ct_by_e (in=in1) dmel.Fb551_si_fusions_unique_flagged;
by fusion_id;
if in1;
rename fusion_id=exonic_region;
run;


proc export data=regression_results
outfile = "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/manuscript/resubmission/Supplemental_Table5_lmm2.csv"
dbms=csv replace;
putnames=yes;
run;

data correlations_of_effects;
merge cegs.corr_ct_by_e(in=in1) dmel.Fb551_si_fusions_unique_flagged;
by fusion_id;
if in1;
rename fusion_id=exonic_region;
run;

proc contents data=correlations_of_effects;
run;


proc export data=correlations_of_effects
outfile = "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/manuscript/resubmission/Supplemental_Table6_lmm2.csv"
dbms=csv replace;
putnames=yes;
run;



