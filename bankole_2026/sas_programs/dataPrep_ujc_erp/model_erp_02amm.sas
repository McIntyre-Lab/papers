
libname lmm_rmg "!MCLAB/lmm_rmg_dros_data/sasdata" ;
libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata" ;
libname sex "!MCLAB/sex_specific_splicing/sasdata" ;

libname mel "!MCLAB/useful_dmel_data/gene_lists/sas_data" ;



/* 


models == ERP_plus

ttest!!


*/


data design_r ;
retain new ;
set lmm_rmg.design_rmg_mel_sim; 
drop sampleID ref techrep;
new = compress("d"||sample);
drop sample ;
rename new = sample ;
new2 = compress("d"||species) ;
drop species;
rename new2 = species ;
run;

proc sort data = design_r nodups ;
by _all_ ;
run;


data design_k ;
set lmm_axk.design_axk_ser; 
drop sampleID ref techrep;
run;

proc sort data = design_k nodups ;
by _all_ ;
run;

%macro model_prep (species) ;

/* for ERPs in genes with readcnts of ge 25 */
data oncall_&species._25 ;
set sex.onCall_gene_&species._on25 ;
if flag_both_sex_off_readCnt25 = 0 ; /* analyze sex-limited and sex-biased */
keep geneID ;
run;

proc sort data = sex.norm_ERP_&species.  ;
by geneID ;
proc sort data = oncall_&species._25;
by geneID ;
run ;

data &species._close  ;
merge oncall_&species._25 (in=in1) sex.norm_ERP_&species. (in=in2) ;
by geneID ;
if in1 ;
run;

data &species._ready ;
set &species._close ;
if log_uq_numTranscripts=. then log_uq_numTranscripts=0;  /*is this right??? */
keep sample ERP_plus log_uq_sumReadNum sex rep geneID ;
run;

%mend ;

%model_prep (dmel) ;
%model_prep (dsim) ;
%model_prep (dser) ;

/*
data dmel_sex ;
set sex.dmel_ujc_cnts_sumtr_sexdet ;
keep dmel6_geneID ;
rename dmel6_geneID = geneID ;
run;

proc sort data = dmel_sex nodups ;
by _all_ ;
proc sort data = dmel_ready ;
by geneID ;
run ;

data dmel_ready_sex ;
merge dmel_ready (in=in1)  dmel_sex (in=in2) ;
by geneID ;
if in2 ;
run;
*/

%macro into (species) ;

proc export data = &species._ready 
outfile = "!MCLAB/sex_specific_splicing/model_output/&species._ready_4_ERP_ttest.csv"
dbms = csv replace ;
run ;
%mend ;

%into (dmel) ;
%into (dsim) ;
%into (dser) ;

%macro ttest (species) ;

proc sort data = &species._ready;
by geneID ERP_plus;
run;

proc ttest data=&species._ready;
by geneID ERP_plus;
class sex;
var log_uq_sumReadNum;
ods output ttests=sex.ttests_sex_&species. equality=sex.eqvr_sex_&species.;
run;

%mend ;

%ttest (dmel) ;
%ttest (dsim) ;
%ttest (dser) ;

/* results file 

geneID
patternID
ttest_equal_var
ttest_unequal_var
pval_equal_var
ttest_var_test
pval_var_test

*/

%macro results (species, genome) ;

proc transpose data = sex.ttests_sex_&species. out = ttest_sex_&species._tvalue prefix=effect_size_;
by geneID ERP_plus ;
id variances ;
var tValue ;
run ;

proc transpose data = sex.ttests_sex_&species. out = ttest_sex_&species._pval prefix=pval_ttest_;
by geneID ERP_plus ;
id variances ;
var probt ;
run ;

proc transpose data = sex.eqvr_sex_&species. out = eqvr_sex_&species._val prefix=probF_;
by geneID ERP_plus ;
id method ;
var Fvalue ;
run ;

proc transpose data = sex.eqvr_sex_&species. out = eqvr_sex_&species._pval prefix=pval_;
by geneID ERP_plus ;
id method ;
var ProbF ;
run ;

proc sort data= ttest_sex_&species._tvalue ;
by geneID ERP_plus ;
proc sort data= ttest_sex_&species._pval ;
by geneID ERP_plus ;
proc sort data= eqvr_sex_&species._val ;
by geneID ERP_plus ;
proc sort data= eqvr_sex_&species._pval  ;
by geneID ERP_plus ;
run ;

data &species._results ;
merge ttest_sex_&species._tvalue ttest_sex_&species._pval eqvr_sex_&species._val eqvr_sex_&species._pval ;
by geneID ERP_plus ;
run ;

data &species._ttest_results ;
retain geneID ERP_plus ERP flagDataOnlyExon flagIR ;
set &species._results ;
ERP = scan(ERP_plus, 1, '|') ;
flagDataOnlyExon = input(scan(ERP_plus, 2, '|'), best32.) ;
flagIR = input(scan(ERP_plus, 3, '|'), best32.)  ;
drop _name_ _label_ ;
run ;

     data WORK.&species._FLAGS    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "/TB20/sex_specific_splicing/fiveSpecies_2_&genome._ujc_er_vs_&species._data_2_&genome._ujc_noMultiGene_flagERP.csv"
     delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat ERP $120. ;
        informat flagDataOnlyExon best32. ;
        informat flagIR best32. ;
        informat geneID $26. ;
        informat seqname $32. ;
        informat strand $1. ;
        informat numJxnHash best32. ;
        informat numAnnotatedER best32. ;
        informat flagNoSkip best32. ;
        informat flagNovel best32. ;
        informat flagERSkip best32. ;
        informat flag5pFragment best32. ;
        informat flag3pFragment best32. ;
        informat flagIntrnlFrgmnt best32. ;
        informat flagFirstER best32. ;
        informat flagLastER best32. ;
        format ERP $120. ;
        format flagDataOnlyExon best12. ;
        format flagIR best12. ;
        format geneID $26. ;
        format seqname $32. ;
        format strand $1. ;
        format numJxnHash best12. ;
        format numAnnotatedER best12. ;
        format flagNoSkip best12. ;
        format flagNovel best12. ;
        format flagERSkip best12. ;
        format flag5pFragment best12. ;
        format flag3pFragment best12. ;
        format flagIntrnlFrgmnt best12. ;
        format flagFirstER best12. ;
        format flagLastER best12. ;
     input
                 ERP  $
                 flagDataOnlyExon
                 flagIR
                 geneID  $
                 seqname  $
                 strand  $
                 numJxnHash
                 numAnnotatedER
                 flagNoSkip
                 flagNovel
                 flagERSkip
                 flag5pFragment
                 flag3pFragment
                 flagIntrnlFrgmnt
                 flagFirstER
                 flagLastER
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

proc sort data = &species._ttest_results ;
by geneID ERP flagDataOnlyExon flagIR ;
proc sort data = &species._flags ;
by geneID ERP flagDataOnlyExon flagIR ;
run;

data &species._ERP_ttest_results1 ;
merge &species._ttest_results (in=in1) &species._flags (in=in2) ;
by  geneID ERP flagDataOnlyExon flagIR ;
if in1 ;
run;

data &species._ERP_ttest_results2 ;
set &species._ERP_ttest_results1 ;
drop numJxnHash numAnnotatedER ;
run;

proc sort data = sex.onCall_erp_&species._flags ;
by geneID ERP_plus ;
proc sort data= &species._ERP_ttest_results2 ; 
by geneID ERP_plus ;
run ;

data &species._ERP_ttest_results ; 
merge &species._ERP_ttest_results2 (in=in1) sex.onCall_erp_&species._flags (in=in2);
by geneID ERP_plus ;
if in1 ;
run;

%mend ;

%results (dmel, dmel6) ;
%results (dsim, dsim2) ;
%results (dser, dser1) ;


/* add dmel6_geneID to sim and ser --> only sex det genes for now */
%macro addmore (species, genome) ;

data link_&species. ;
set sex.fivespecies_&genome._w_geneID_02amm ;
keep &genome._geneID dmel6_geneID ;
run;

proc sort data=  link_&species.  nodups;
by _all_ ;
run;

data &species._ERP_ttest_results_A ;
set &species._ERP_ttest_results ;
rename geneID = &genome._geneID ;
run ;

data &species._ERP_ttest_results_B ;
merge link_&species. (in=in1) &species._ERP_ttest_results_A (in=in2);
by &genome._geneID ;
if in2 ;
run ;

%mend ;

%addmore (dsim, dsim2)  ;
%addmore (dser, dser1)  ;

data dmel_ERP_ttest_results_B ;
set dmel_ERP_ttest_results ;
rename geneID = dmel6_geneID ;
run;

/* merge in gene symbols */
proc import datafile = "!MCLAB/useful_dmel_data/flybase650/dmel_annotation/fbgn_annotation_ID.csv"
out = dmel_symbols 
dbms = csv replace ;
guessingrows = MAX ;
run ;

data dmel_symbols ;
set  dmel_symbols ;
keep symbol primary_FBgn ;
rename primary_FBgn = dmel6_geneID ;
run;

proc sort data= dmel_symbols nodups ;
by dmel6_geneID symbol ;
run;

data sex.dmel_ERP_ttest_results ;
merge dmel_symbols (in=in1) dmel_ERP_ttest_results_B (in=in2) ;
by dmel6_geneID ;
if in2 ;
run ;

proc export data = sex.dmel_ERP_ttest_results 
outfile = "!MCLAB/sex_specific_splicing/model_output/dmel_ERP_ttest_results.csv"
dbms = csv replace ;
run ;

/* add sim gene symbols */
libname sim "!MCLAB/useful_dsim_data/flybase202/sas_data";

data dsim_symbols ; 
set sim.dsim_genes2go_ortho ;
keep symbol FBgn ;
rename FBgn = dsim2_geneID ;
run;

proc sort data= dsim_symbols nodups ;
by dsim2_geneID symbol ;
run;

data sex.dsim_ERP_ttest_results ;
merge dsim_symbols (in=in1) dsim_ERP_ttest_results_B (in=in2) ;
by dsim2_geneID ;
if in2 ;
run ;

proc export data = sex.dsim_ERP_ttest_results 
outfile = "!MCLAB/sex_specific_splicing/model_output/dsim_ERP_ttest_results.csv"
dbms = csv replace ;
run ;

/* no link for dser symbols */
data sex.dser_ERP_ttest_results  ;
set dser_ERP_ttest_results_B  ;
run ;

proc export data = sex.dser_ERP_ttest_results 
outfile = "!MCLAB/sex_specific_splicing/model_output/dser_ERP_ttest_results.csv"
dbms = csv replace ;
run ;



