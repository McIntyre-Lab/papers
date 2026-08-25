
PROC IMPORT OUT= sim_gene 
           DATAFILE="~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/summary_files/gene_level_summary_from_dataFile_jxnHash_w_genesetid_dsim.csv" 
      DBMS=CSV REPLACE; 
     GETNAMES=YES; 
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN; 


PROC IMPORT OUT= sim_ujc 
           DATAFILE="~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/datafiles/datafile_jxnHash_dsim2.csv" 
      DBMS=CSV REPLACE; 
     GETNAMES=YES; 
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN; 

data sim_both;
set sim_gene;
if sexBias="B" ;
num_ujc_F=num_ujc_limited_F+num_ujc_bias_F;
num_ujc_M=num_ujc_limited_M+num_ujc_bias_M;
flag_ujc_bias=0;
if num_ujc_F>0 and num_ujc_M>0 then flag_ujc_bias=1;

run;

proc freq data=sim_both;
tables flag_ujc_bias;
run;

/* of the  608  genes*/

/*note some genes are bias in both because of ERP not UJC we do not expect all genes */

data both_w_ujc;
set sim_both;
if flag_ujc_bias=1;
keep geneID;
run;

/*435*/
proc sort data=both_w_ujc;
by geneID;

proc sort data=sim_ujc;
by geneID;

/*only significant ujc*/


data find_ujc_inboth;
merge sim_ujc (in=in1) both_w_ujc(in=in2);
by geneID;
if in2;
if sexClass="unanalyzed" or sexClass="unbiased" then delete;
run;
/*ujc have sexbias in these genes of these  in this list*/


data female_bias_ujc_inboth;
set find_ujc_inboth;
if flag_sex_bias_F=1 or  flag_sex_limited_F_rc50=1;
run;
/*1158*/
data male_bias_ujc_inboth;
set find_ujc_inboth;
if flag_sex_bias_M=1 or  flag_sex_limited_M_rc50=1;
run;
/*874*/
proc freq data=female_bias_ujc_inboth;
tables structural_category*flagnovel_erp;
run;


proc freq data=male_bias_ujc_inboth;
tables structural_category*flagnovel_erp;
run;

/*very few novel ERP-  in male and  in female*/
/*NIC in female,  in male  fsm in male and fsm in female*/

/*check against gene_summary*/

proc freq data=female_bias_ujc_inboth;
tables geneID/out=count_sim_female;

data count_sim_female1;
set count_sim_female;
rename count=sim_num_ujc_f;
keep geneID count;

proc freq data=male_bias_ujc_inboth;
tables geneID/out=count_sim_male;
run;

data count_sim_male1;
set count_sim_male;
rename count=sim_num_ujc_m;
keep geneID count;


proc sort data=sim_both;
by geneID;

proc sort data=count_sim_female1;
by geneID;

proc sort data=count_sim_male1;
by geneID;

data check_counts oops;
merge sim_both(in=in1) count_sim_female1(in=in2) count_sim_male1(in=in3);
by geneID;
if in1 and in2 and in3 then output check_counts;
 else output oops;
run;
/*435 genes */

data check_numujc;

set check_counts;
if sim_num_ujc_m = (num_ujc_bias_M+num_ujc_limited_M) then flag_male_match=1;
	else flag_male_match=0;
if sim_num_ujc_f = (num_ujc_bias_F+num_ujc_limited_F) then flag_female_match=1;
	else flag_female_match=0;

proc freq data=check_numujc;
table flag_male_match flag_female_match;
run;
/*matches!

/*this should have zero oops*/
data check_counts oops;
merge count_sim_female1(in=in2) count_sim_male1(in=in3);
by geneID;
if  in2 and in3 then output check_counts;
 else output oops;
run;
/*yeah!*/

PROC EXPORT DATA= male_bias_ujc_inboth
            OUTFILE= "!MCLAB/sex_specific_splicing/TranD_sexBias/dsim_male_ujc_4_TranD.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


PROC EXPORT DATA= female_bias_ujc_inboth
            OUTFILE= "!MCLAB/sex_specific_splicing/TranD_sexBias/dsim_female_ujc_4_TranD.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




