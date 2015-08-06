libname sem 'S:\McIntyre_Lab\cegs_sem_sd_paper\sas_data';
libname SEM '!MCLAB/cegs_sem_sd_paper/sas_data';

data look_flags;
set sem.validation_set;
run;

proc freq data=look_flags;
tables flag_sem_added_gene*(flag_dsxNull_repressed
flag_chang_ds_tra
flag_chang_ups_tra
flag_tra_bs
flag_goldman_ds_tra);
run;

data check;
set look_flags;
if flag_tra_bs=1 and flag_sem_added_gene=1;
run;

proc freq data=check;
tables flag_sem_added_gene*(flag_dsxNull_repressed
flag_chang_ds_tra
flag_chang_ups_tra
flag_tra_bs
flag_goldman_ds_tra);
run;

*merged dataset from cegsV_export_adding_genes_nofilter_yp2.sas;

%macro pull_models(fbgn,gene);

data &gene;
set merged;
where gene ?"&gene" or gene ?"&fbgn";
run;

proc sort data=&gene;
by bic;
run;
%mend;
%pull_models (FBgn0261243,psa);
%pull_models (FBgn0005616,msl2);
%pull_models (FBgn0259481,mob); *yes;
%pull_models (FBgn0004587,B52); *yes;
%pull_models (FBgn0263396,sqd); *yes;
%pull_models (FBgn0014870,psi);*no;
%pull_models (FBgn0262737,Mub); *not in merged?;
%pull_models (FBgn0003261,Rm62); *winner;
%pull_models (FBgn0261552,ps); *no evidence;
%pull_models (FBgn0263391,hts); *downstream sxl;
%pull_models (FBgn0002567,ltd); *wahoo;
%pull_models (FBgn0264560,garz); *no;
%pull_models (FBgn0027654,jdp); *no;
%pull_models (FBgn0033844,CG6016); *no;
%pull_models (FBgn0038735,CG4662); *no;

proc freq data=look_flags;
tables flag_sem_added_gene*flag_chang_ds_tra*flag_chang_ups_tra;
run;

proc freq data=look_flags;
where flag_sem_added_gene=1 and (flag_chang_ds_tra=1 or flag_chang_ups_tra=1);
tables flag_tra_bs;
run;


proc freq data=look_flags;
where flag_sem_added_gene=1 and (flag_chang_ds_tra=1 or flag_chang_ups_tra=1);
tables flag_goldman_ds_tra;
run;


data check2;
set look_flags;
if flag_sem_added_gene=1 and flag_goldman_ds_tra=1 and (flag_chang_ds_tra=1 or flag_chang_ups_tra=1);
run;


*from useful dmel gene lists import mcintyre 2006;

data splice;
set lookup_fbgn;
if psex_by_probe ne .;
run;

data list_splice;
set splice;
rename fbgn=gene;
keep fbgn;

proc sort data=list_splice nodups out=fbgn_splice_nodup;
by gene;
run;

proc sort data=list_splice;
by gene;
proc sort data=merged;
by gene;

data check_model_4_splice;
merge merged (in=in1) list_splice(in=in2);
by gene;
if in2;
run;
*sloppy merge;

data no_match match;
set check_model_4_splice;
if model='' then output no_match;
else output match;
run;

proc sort data=match;
by gene bic;
run;

proc freq data=match;
tables gene/out=count_gene;
run;

data check_best;
set match;
by gene;
if first.gene;
run;

proc freq data=check_best;
tables model;
run;





/*
Fear 2014:
flag_sem_added_gene = 1 if the SEM model added the gene to the SD pathway
flag_dsxNull_induced = 1 if a gene had increased expression when dsx was knocked out
flag_dsxNull_repressed= 1 if a gene had decreased expression when dsx was knocked out (since dsx is a TF I could guess these are more likely direct targets)
 
Chang 2011:
flag_chang_female_bias = 1 if there was a diff between MvsF and it was biased toward females
flag_chang_male_bias= 1 if there was a diff between MvsF and it was biased toward males
flag_chang_tra_bias = 1 if there was a diff between FvsTra pseudo males and it was biased toward tra pseudomales
flag_chang_traf_bias = 1 if there was a diff between FvsTra pseudo males and it was biased toward females
flag_chang_ds_tra = 1 if looking at the results the gene appears to be downstream of tra
flag_chang_ups_tra = 1 if looking at the results the gene appears to be upstream of tra
flag_tra_bs = 1 if contains a tra dna binding site
 
Goldman 2007:
flag_goldman_female_bias = 1 if there was a diff between MvsF and it was biased toward females
flag_goldman_male_bias bias = 1 if there was a diff between MvsF and it was biased toward males
flag_goldman_ds_tra =  1 if looking at the results the gene appears to be downstream of tra
*/
