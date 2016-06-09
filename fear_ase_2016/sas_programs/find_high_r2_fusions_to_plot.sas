proc datasets library=work kill;
run; quit;

      data WORK.var_stack    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile 'Y:\cegs_ase_paper\output\cis_trans_variance_data_stack.csv' delimiter = ','
 MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat type $5. ;
         informat fusion_id $9. ;
         informat mating_status $1. ;
         informat variance best32. ;
         format type $5. ;
         format fusion_id $9. ;
         format mating_status $1. ;
         format variance best12. ;
      input
                 type $
                  fusion_id $
                 mating_status $
                 variance
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

      data WORK.r2_stack    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile 'Z:\cegs_ase_paper\output\cis_trans_r2_data_stack.csv' delimiter = ',' MISSOVER
 DSD lrecl=32767 firstobs=2 ;
        informat type $15. ;
        informat fusion_id $9. ;
        informat mating_status $1. ;
        informat R2 best32. ;
        format type $15. ;
        format fusion_id $9. ;
         format mating_status $1. ;
         format R2 best12. ;
     input
                 type $
                 fusion_id $
                 mating_status $
                 R2
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

/* For R2_cis */
data r2_cis;
set r2_stack;
if type="R2_cis";
rename variance= cis_variance;
drop type;
run;

data var_cis;
set var_stack;
if type="cis";
rename R2= R2_cis;
drop type;
run;

proc sort data=r2_cis;
by fusion_id mating_status;
proc sort var_cis;
by fusion_id mating_status;
run;

data r2_cis_vs_var;
merge r2_cis var_cis;
by fusion_id mating_status;
run;

data high;
set r2_cis_vs_var;
if R2 < 0.95 then delete;
run;

proc sort data=high;
by variance;
run;

/* For R2_trans */
data r2_trans;
set r2_stack;
if type="R2_diff_trans";
drop type;
run;

data var_trans;
set var_stack;
if type="trans";
drop type;
run;

proc sort data=r2_trans;
by fusion_id mating_status;
proc sort data=var_trans;
by fusion_id mating_status;
run;

data r2_trans_vs_var;
merge r2_trans var_trans;
by fusion_id mating_status;
run;

data high_trans;
set r2_trans_vs_var;
if R2 < 0.80 then delete;
run;

proc sort data=high_trans;
by variance;
run;

/* For int */
data r2_int;
set r2_stack;
if type="R2_diff_int";
drop type;
run;

proc sort data=r2_int;
by fusion_id mating_status;
proc sort data=var_cis;
by fusion_id mating_status;
run;

data r2_int_vs_var;
merge r2_int var_cis;
by fusion_id mating_status;
run;

data high_int;
set r2_int_vs_var;
if R2 < 0.50 then delete;
run;

proc sort data=high_int;
by variance;
run;


/* Gene lists for high R2 cis, R2 trans, R2 interaction */
proc sort data=r2_cis;
by descending R2;
proc sort data=r2_trans;
by descending R2;
proc sort data=r2_int;
by descending R2;
run; *1760 obs in each dataset;

data r2_cis_top;
set r2_cis ;
if R2 > 0.9 ;
flag_cis = 1;
rename r2=r2_cis;
run;

data r2_trans_top;
set r2_trans;
if R2 > 0.4 ;
flag_trans=1;
rename r2=r2_trans;
run;

data r2_int_top;
set r2_int;
if R2 > 0.4;
flag_int=1;
rename r2=r2_int;
run;

proc sort data=r2_cis_top;
by fusion_id mating_status;
proc sort data=r2_trans_top;
by fusion_id mating_status;
proc sort data=r2_int_top;
by fusion_id mating_status;
run;

data overlap unique;
merge r2_cis_top (in=in1) r2_trans_top (in=in2) r2_int_top (in=in3);
by fusion_id mating_status;
if in1 and in2 and in3 then output overlap;
else output unique;
run; *0 obs in overlap, 435 obs in unique;

proc sort data=r2_cis_top out=nodup_r2_cis nodupkey;
by fusion_id;
proc sort data=r2_trans_top out=nodup_r2_trans nodupkey;
by fusion_id;
proc sort data=r2_int_top out=nodup_r2_int nodupkey;
by fusion_id;
run;

*separate and flag if M or V;
%macro sep (type);
data r2_cis_m;
set r2_cis_top;
if mating_status="M";
flag_M = 1;
drop mating_status;
run;

data r2_cis_V;
set r2_cis_top;
if mating_status="V";
flag_V = 1;
drop mating_status;
run;
%mend;
%sep (int);
%sep (trans);
%sep (cis);


%macro sort (type);
proc sort data=r2_&type._V;
by fusion_id ;
proc sort data=r2_&type._M;
by fusion_id ;
run;

data &type._sbs;
merge r2_&type._M r2_&type._V;
by fusion_id ;
if flag_V ne 1 then flag_V = 0;
if flag_M ne 1 then flag_M = 0;
flag_both = flag_M + flag_V ;
run;

proc sort data=&type._sbs;
by fusion_id;
run;

%mend;
%sort (int);
%sort (trans);
%sort (cis);

libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
proc sort data=dmel.fb551_si_fusions_unique;
by fusion_id;
run;

data cis_w_gene;
merge cis_sbs (in=in1) dmel.fb551_si_fusions_unique;
by fusion_id;
if in1;
run;

data trans_w_gene;
merge trans_sbs (in=in1) dmel.fb551_si_fusions_unique;
by fusion_id;
if in1;
run;

data int_w_gene;
merge int_sbs (in=in1) dmel.fb551_si_fusions_unique;
by fusion_id;
if in1;
run;

/*output the datasets*/
data output_cis;
set cis_w_gene;
keep fusion_id symbol_cat flag_M flag_V flag_both;
run;

data output_trans;
set trans_w_gene;
keep fusion_id symbol_cat flag_M flag_V flag_both;
run;

data output_int;
set int_w_gene;
keep fusion_id symbol_cat flag_M flag_V flag_both;
run;

proc sort data=output_int;
by symbol_cat;
proc sort data=output_trans;
by symbol_cat;
proc sort data=output_cis;
by symbol_cat;
run;

*get unique only;
data oops overlap1 overlap2 overlap3 int_only trans_only cis_only;
merge output_int (in=in1) output_trans (in=in2) output_cis (in=in3);
by symbol_cat;
if in1 and in2 then output overlap1; *4 overlap int and trans;
else if in1 and in3 then output overlap2; *3 overlap int and cis;
else if in2 and in3 then output overlap3; *12 overlap cis and trans;
else if in1 then output int_only; *22 obs;
else if in2 then output trans_only; *51 obs;
else if in3 then output cis_only; *175 obs;
else output oops;
run; *0 in overlap;

proc export data=cis_only
outfile = "Z:/cegs_ase_paper/output/gene_categories/top_r2_cis_genes.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=trans_only
outfile = "Z:/cegs_ase_paper/output/gene_categories/top_r2_trans_genes.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=int_only
outfile = "Z:/cegs_ase_paper/output/gene_categories/top_r2_interaction_genes.csv"
dbms=csv replace;
putnames=yes;
run;

* find R2 values for graphs ;
data men_R2;
set r2_stack;
if fusion_id = "S52856_SI" and mating_status = "V";
run;

data chop24_R2;
set r2_stack;
if fusion_id = "S13419_SI" and mating_status = "V";
run;

data CG3036_R2;
set r2_stack;
if fusion_id = "S2713_SI" and mating_status = "V";
run;
