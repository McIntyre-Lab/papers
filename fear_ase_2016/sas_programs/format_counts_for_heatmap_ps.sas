libname cegs "McLab/cegs_ase_paper/sas_data";
filename mymacros 'McLab/maize_ozone/2014/sas_analysis/macros';
options SASAUTOS=(sasautos mymacros);

       data WORK.design_file    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
       infile 'McLab/cegs_ase_paper/output/CEGS_49_lines.txt' delimiter='09'x MISSOVER DSD
 lrecl=32767 firstobs=1 ;
          informat col1 $7. ;
          format col1 $7. ;
       input
                   col1 $
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;

*make new columns based on length of line name ;
data new_design1;
set design_file (firstobs=1 obs=21);
line= substr(put(col1,10.), 1,4);
tech_rep = substr(put(col1,10.), 6,2);
drop col1;
run;

data new_design2;
set design_file (firstobs=22 obs=27);
line= substr(put(col1,10.), 1,3);
tech_rep = substr(put(col1,10.), 5,2);
drop col1;
run;

data new_design3;
set design_file (firstobs=28 obs=214);
line= substr(put(col1,10.), 1,4);
tech_rep = substr(put(col1,10.), 6,2);
drop col1;
run;

data new_design4;
set design_file (firstobs=215 obs=220);
line= substr(put(col1,10.), 1,3);
tech_rep = substr(put(col1,10.), 5,2);
drop col1;
run;

data new_design5;
set design_file (firstobs=221 obs=232);
line= substr(put(col1,10.), 1,4);
tech_rep = substr(put(col1,10.), 6,2);
drop col1;
run;

data new_design6;
set design_file (firstobs=233 obs=238);
line= substr(put(col1,10.), 1,3);
tech_rep = substr(put(col1,10.), 5,2);
drop col1;
run;

data new_design7;
set design_file (firstobs=239 obs=250);
line= substr(put(col1,10.), 1,4);
tech_rep = substr(put(col1,10.), 6,2);
drop col1;
run;

data new_design8;
set design_file (firstobs=251 obs=310);
line= substr(put(col1,10.), 1,3);
tech_rep = substr(put(col1,10.), 5,2);
drop col1;
run;

*set all design files together;
data design_file2;
set new_design: ;
run; *310 obs;

*import each dataset;
%macro import (geno, rep);
       data &geno._&rep    ;
       %let _EFIERR_ = 0; 
       infile "McLab/cegs_ase_paper/output/counts_gene_ps/&geno._&rep..txt" delimiter='09'x MISSOVER
 DSD lrecl=32767 firstobs=1 ;
          informat col1 $30. ;
          format col1 $30. ;
       input
                   col1 $
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  
	       run;

data to_stack_&rep._&geno ;
set &geno._&rep;
chrom= substr(put(col1,10.), 1,2);
pos = substr(put(col1,20.), 4,7);
count_1 = substr(put(col1,20.), 12,2);
count = count_1 + 0;
genotype = "&geno";
rep = "&rep";
drop col1;
run;

%mend;
%iterdataset(dataset=design_file2, function=%nrstr(%import(&line, &tech_rep);));

*drop all empty datasets-- did not have bam files for these;
proc sql noprint ;
select 'drop table work.'||memname into :empties separated by ';'
from dictionary.tables
where libname='WORK' and nobs=0
;
&empties;
quit;

*stack all count files together then proc freq to see which genotype*reps need added together;
data mated_count_stack;
set to_stack_m: ;
run;

data virgin_count_stack;
set to_stack_v: ;
run;

*proc freq to update design file;
proc sort data=mated_count_stack out=mated_nodup nodupkey;
by genotype rep;
run;
proc freq data=mated_nodup;
tables genotype*rep / out=mated_design_file;
run;

proc sort data=virgin_count_stack out=virgin_nodup nodupkey;
by genotype rep;
run;
proc freq data=virgin_nodup;
tables genotype*rep / out=virgin_design_file;
run;

proc datasets noprint;
delete to_stack_: ;
delete r: ;
delete w: ;
run; quit;

/*means for each position, chrom, genotype*/
proc means data=virgin_count_stack ;
var count;
class chrom pos genotype;
output out=virgin_mean mean=mean_count sum=sum_count max=max_count;
run;

data virgin_max;
set virgin_mean;
if _TYPE_ = 1;
keep genotype max_count;
run;

data virgin_mean2;
set virgin_mean;
if _type_ = 7;
drop _TYPE_ _FREQ_ chrom sum_count max_count;
run;

proc sort data=virgin_mean2;
by genotype;
proc sort data=virgin_max;
by genotype;
run;

data virgin_merge;
merge virgin_max virgin_mean2;
by genotype;
run;

*scaled by the maximum count value per line;
data virgin_scaled;
set virgin_merge;
scaled_count = mean_count / max_count;
drop mean_count max_count;
run;

/*MATED*/
proc means data=mated_count_stack ;
var count;
class chrom pos genotype;
output out=mated_mean mean=mean_count sum=sum_count max=max_count;
run;

data mated_mean2;
set mated_mean;
if _type_ = 7;
drop _TYPE_ _FREQ_ chrom max_count sum_count;
run;

data mated_max;
set mated_mean;
if _TYPE_ = 1;
keep genotype max_count;
run;

data mated_mean2;
set mated_mean;
if _type_ = 7;
drop _TYPE_ _FREQ_ chrom sum_count max_count;
run;

proc sort data=mated_mean2;
by genotype;
proc sort data=mated_max;
by genotype;
run;

data mated_merge;
merge mated_max mated_mean2;
by genotype;
run;

*scaled by the maximum count value per line;
data mated_scaled;
set mated_merge;
scaled_count = mean_count / max_count;
drop mean_count max_count;
run;

*export the long form and try to manipulate in R;
proc export data=mated_mean2
outfile = "McLab/cegs_ase_paper/output/mated_ps_counts_for_heatmap.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=virgin_mean2
outfile = "McLab/cegs_ase_paper/output/virgin_ps_counts_for_heatmap.csv"
dbms=csv replace;
putnames=yes;
run;

*export scaled data;
proc export data=mated_scaled
outfile = "McLab/cegs_ase_paper/output/mated_ps_counts_scaled.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=virgin_scaled
outfile = "McLab/cegs_ase_paper/output/virgin_ps_counts_scaled.csv"
dbms=csv replace;
putnames=yes;
run;

proc freq data=virgin_scaled;
tables scaled_count / out=chk_scale;
run;

/************* STOP *************/

/* TRANSPOSE IN SAS IS NOT GETTING ALL OF THE OBSERVATIONS!!!! */
*we need positions across the top;
proc transpose data=mated_mean2 out=mated_wide prefix=ch3R_ ;
by genotype;
id pos;
var mean_count;
run;

data mated_wide;
set mated_wide;
drop _NAME_;
run;

*change all . to 0-- no coverage;
data mated_wide2;
set mated_wide;
array a(*) _numeric_;
do i=1 to dim(a);
if a(i) = . then a(i) = 0;
end;
drop i;
run;

*virgin;
proc transpose data=virgin_mean2 out=virgin_wide prefix=ch3R_ ;
by genotype;
id pos;
var mean_count;
run;

data virgin_wide;
set virgin_wide;
drop _NAME_;
run;

*change all . to 0-- no coverage;
data virgin_wide2;
set virgin_wide;
array a(*) _numeric_;
do i=1 to dim(a);
if a(i) = . then a(i) = 0;
end;
drop i;
run;

*proc export doesn't export the entire file-- use this command-- it takes FOREVER;
ods csv file = "McLab/cegs_ase_paper/output/virgin_ps_count2.csv";
proc report data=virgin_wide2 nowd;
define _all_ / display;
run;
ods csv close;

ods csv file = "McLab/cegs_ase_paper/output/mated_ps_count2.csv";
proc report data=mated_wide2 nowd;
define _all_ / display;
run;
ods csv close;


