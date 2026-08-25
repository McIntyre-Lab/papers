
libname lmm_rmg "!MCLAB/lmm_rmg_dros_data/sasdata" ;
libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata" ;
libname sex "!MCLAB/sex_specific_splicing/sasdata" ;

libname mel "!MCLAB/useful_dmel_data/gene_lists/sas_data" ;



/* 


models == UJC

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

/* for UJC in genes with readcnts of ge 25 */
data oncall_&species._25 ;
set sex.onCall_gene_&species._on25 ;
if flag_both_sex_off_readCnt25 = 0 ; /* analyze sex-limited and sex-biased */
keep geneID ;
run;

proc sort data = sex.norm_UJC_&species.  ;
by geneID ;
proc sort data = oncall_&species._25 ;
by geneID ;
run ;

data &species._close  ;
merge oncall_&species._25 (in=in1) sex.norm_UJC_&species. (in=in2) ;
by geneID ;
if in1 ;
run;

data &species._ready ;
set &species._close ;
if log_uq_numTranscripts=. then log_uq_numTranscripts=0;  /*is this right??? */
keep sample jxnHash log_uq_sumReadNum sex rep geneID ;
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

/*
%macro into (species) ;

proc export data = &species._ready 
outfile = "!MCLAB/sex_specific_splicing/model_output/&species._ready_4_ERP_ttest.csv"
dbms = csv replace ;
run ;
%mend ;

%into (dmel) ;
%into (dsim) ;
%into (dser) ;
*/
%macro ttest (species) ;

proc sort data = &species._ready;
by jxnHash;
run;

proc ttest data=&species._ready;
by jxnHash;
class sex;
var log_uq_sumReadNum;
ods output ttests=sex.ttests_sex_ujc_&species. equality=sex.eqvr_sex_ujc_&species.;
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

 
/* combine ujc ttest results

input:
    work.&species._ttest_results from model_ujc
    sex.cross_species_&var._2_dmel6


*/




