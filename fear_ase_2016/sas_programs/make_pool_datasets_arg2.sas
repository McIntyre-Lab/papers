libname cegs "Z:/cegs_ase_paper/sas_data";

PROC IMPORT OUT= WORK.hubs 
            DATAFILE= "Z:/cegs_ase_paper/Pool_2015_gene_lists/Table7_gene_list.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	guessingrows=1689;
RUN;

*i get 132 from the number of gene columns;
data gene_key;
set hubs (firstobs=1 obs=132);
number_key= _n_ ;
keep number_key;
run;

data drop_hubs;
set hubs;
row= _n_;
rename Genes_in_outlier_regions= gene132;
drop GO_number flag_differentiation_hub ;
run;

filename mymacros 'Z:/maize_ozone/2014/sas_analysis/macros';
options SASAUTOS=mymacros;
%include 'Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas';


data hubs2;
set hubs;
rename Genes_in_outlier_regions= gene0;
run;

%macro by_gene (ID);

data gene_&ID;
set hubs2;
keep gene&ID;
run;

data clean_gene_&ID;
set gene_&ID;
if gene&ID = '' then delete;
rename gene&ID = table7_gene;
run;

%mend;
%iterdataset(dataset=gene_key, function=%nrstr(%by_gene(&number_key);));

data table7_genes;
set clean_gene_:;
flag_table7_gene=1;
run;

proc sort data=table7_genes nodupkey;
by table7_gene;
run; *140 genes;

proc datasets noprint;
delete clean_gene_: ;
delete gene: ;
run; quit;

/* TABLE 4 EUROPEAN GENES */
PROC IMPORT OUT= WORK.european 
            DATAFILE= "Z:/cegs_ase_paper/Pool_2015_gene_lists/Table4_european_gene_list.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	guessingrows=3389;
RUN;

*i get 132 from the number of gene columns;
data gene_key;
set european (firstobs=1 obs=757);
number_key= _n_ ;
keep number_key;
run;

*rename gene0 since count starts at 1 ;
data european2;
set european;
rename gene0= gene757;
run;

%macro by_gene (ID);

data gene_&ID;
set european2;
keep gene&ID;
run;

data clean_gene_&ID;
set gene_&ID;
if gene&ID = '' then delete;
table4_european_gene =put(gene&ID, 20.);
run;

%mend;
%iterdataset(dataset=gene_key, function=%nrstr(%by_gene(&number_key);));

data table4_european_genes;
set clean_gene_:;
flag_table4_europe_gene=1;
keep flag_table4_europe_gene table4_european_gene;
run;

proc sort data=table4_european_genes nodupkey;
by table4_european_gene;
run; *850 obs;

proc datasets noprint;
delete clean_gene_: ;
delete gene: ;
run; quit;


/* TABLE 4 AFRICAN GENES */
PROC IMPORT OUT= WORK.african 
            DATAFILE= "Z:/cegs_ase_paper/Pool_2015_gene_lists/Table4_african_gene_list.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	guessingrows=3344;
RUN;

*i get 568 from the number of gene columns;
data gene_key;
set african (firstobs=1 obs=568);
number_key= _n_ ;
keep number_key;
run;

*rename gene0 since count starts at 1 ;
data african2;
set african;
rename gene0= gene568;
run;

%macro by_gene (ID);

data gene_&ID;
set african2;
keep gene&ID;
run;

data clean_gene_&ID;
set gene_&ID;
if gene&ID = '' then delete;
table4_african_gene =put(gene&ID, 20.);
run;

%mend;
%iterdataset(dataset=gene_key, function=%nrstr(%by_gene(&number_key);));

data table4_african_genes;
set clean_gene_:;
flag_table4_african_gene=1;
keep flag_table4_african_gene table4_african_gene;
run;

proc sort data=table4_african_genes nodupkey;
by table4_african_gene;
run; *639 genes;

proc datasets noprint;
delete clean_gene_: ;
delete gene: ;
run; quit;

data table4_african_genes ;
set table4_african_genes ;
rename table4_african_gene= gene_name;
run;

data table4_african_genes1 ;
set table4_african_genes (firstobs=2 obs=639);
run;

data table4_european_genes;
set table4_european_genes;
rename table4_european_gene=gene_name;
run; *850 obs;

data table7_genes;
set table7_genes;
rename table7_gene= gene_name;
run;

proc sort data= table4_african_genes1;
by gene_name;
proc sort data=table4_european_genes;
by gene_name;
proc sort data=table7_genes;
by gene_name;
run;

data almost;
length gene_name $ 25;
merge table4_african_genes1 table4_european_genes table7_genes;
by gene_name;
if flag_table4_african_gene = . then flag_table4_african_gene=0;
if flag_table4_europe_gene= . then flag_table4_europe_gene=0;
if flag_table7_gene = . then flag_table7_gene=0;
run; *1594 obs;

*one gene is not behaving;
data extra_gene;
length gene_name $ 25;
gene_name="15219";
flag_table4_african_gene=1;
flag_table4_europe_gene=0;
flag_table7_gene=0;
run;

data cegs.pool_2015_gene_lists ;
length gene_name $ 25;
set almost extra_gene;
run;


/****** RITA DATASETS *******/
libname rita "McLab/berlin-c167-hyb-solexa/Analysis RMG/SAS Data/Bayesian";

data defense;
set rita.gene_results_2_def;
keep chrom gene_symbol sim_biased mel_biased bias_dir_sim1 Sackton_defense defense_gene our_defense;
run;

data sex;
set rita.gene_results_2_sex;
keep chrom gene_symbol sim_biased mel_biased bias_dir_sim1 sex_differential_gene downstream_of_tra sex_bias_opposite_tra dsx_set regulated_tra_not_dsx_not_frum regulated_only_one_dsx fru_head_set fru_CNS_head_set sexreg_gene fru_set tra_set alldsx_set;
run;

data all_by_gene;
set rita.all_results_bygene_and_mktests;
keep gene_symbol chrom  sim_biased mel_biased bias_dir_sim1 inter5primeofgene_adptv inter3primeofgene_adptv MKsig_any MK_adaptive_any inter5primeofgene_sig inter3primeofgene_sig ;
run; 

data gene_results_keyword;
set rita.gene_results_2_go_2010_keywords;
keep gene_symbol chrom sim_biased mel_biased bias_dir_sim1 vision_gene nervous_system_gene olfaction_gene silencing_gene germline_gene;
run;

proc sort data=defense;
by gene_symbol;
proc sort data=sex;
by gene_symbol;
proc sort data=gene_results_keyword;
by gene_symbol;
proc sort data=all_by_gene;
by gene_symbol;
run;

data cegs.rita_results_w_flags;
merge defense sex gene_results_keyword all_by_gene;
by gene_symbol;
run; *7515 obs;

proc export data=cegs.rita_results_w_flags
outfile= "McLab/cegs_ase_paper/rita_gene_lists.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=cegs.pool_2015_gene_lists
outfile= "McLab/cegs_ase_paper/pool_2015_gene_lists.csv"
dbms=csv replace;
putnames=yes;
run;
