/*

all ERP data files

*/


libname melsas "!MCLAB/useful_dmel_data/flybase650/sas_data" ;

libname sex "!MCLAB/sex_specific_splicing/sasdata" ;


/* wide format file with raw erp counts */
%macro erp (species) ;
proc sort data = sex.&species._erp_cnts_sumtr_stack_w0s ;
by erp_plus geneID;
run;

proc transpose data = sex.&species._erp_cnts_sumtr_stack_w0s  out = &species._raw_erp_wide prefix = erpCnts_;
by erp_plus geneID;
id sample ;
var erp_sumReadNum ;
run ;

data &species._erp_wide2 ;
set &species._raw_erp_wide  ;
drop _name_ ;
run;
%mend ;

%erp (dmel) ;
%erp (dsim) ;
%erp (dser) ;
%erp (dyak) ;
%erp (dsan) ;

/* create wide format file with log_uq_sumReadNum for merging in  */
%macro uqs (species, genome) ;

proc sort data= sex.norm_erp_&species. ;
by erp_plus geneID;
run;

proc transpose data = sex.norm_erp_&species. out = &species._uq_erp_wide prefix=logUQ_;
by erp_plus geneID;
id sample ;
var log_uq_sumReadNum ;
run;

data &species._uq_erp_wide2 ;
set &species._uq_erp_wide;
drop _name_ ;
run ;
%mend ;

%uqs (dmel, dmel6) ;
%uqs (dsim, dsim2) ;
%uqs (dser, dser1) ;


%macro onoffs (species ) ;
/* add on off flags - includes flag analyzable */
proc sort data = sex.onCall_erp_&species._flags ;
by erp_plus geneID ;
run ;

%mend ;

%onoffs (dmel) ;
%onoffs (dsim) ;
%onoffs (dser) ;
%onoffs (dyak) ;
%onoffs (dsan) ;


/* ttest results */
%macro ttest (species) ;

proc sort data = sex.ttests_sex_erp_&species. ;
by erp_plus geneID ;
proc sort data = sex.eqvr_sex_erp_&species. ;
by erp_plus geneID ;
run;

proc transpose data = sex.ttests_sex_erp_&species. out = ttests_sex_erp_&species._tvalue prefix=effect_size_;
by erp_plus geneID ;
id variances ;
var tValue ;
run ;

proc transpose data = sex.ttests_sex_erp_&species. out = ttests_sex_erp_&species._pval prefix=pval_ttest_;
by erp_plus geneID ;
id variances ;
var probt ;
run ;

proc transpose data = sex.eqvr_sex_erp_&species. out = eqvr_sex_erp_&species._val prefix=probF_;
by erp_plus geneID ;
id method ;
var Fvalue ;
run ;

proc transpose data = sex.eqvr_sex_erp_&species. out = eqvr_sex_erp_&species._pval prefix=pval_;
by erp_plus geneID ;
id method ;
var ProbF ;
run ;

proc sort data= ttests_sex_erp_&species._tvalue ;
by erp_plus geneID ;
proc sort data= ttests_sex_erp_&species._pval ;
by erp_plus geneID ;
proc sort data= eqvr_sex_erp_&species._val ;
by erp_plus geneID ;
proc sort data= eqvr_sex_erp_&species._pval ;
by erp_plus geneID ;
run ;

data &species._results_erp ;
merge   ttests_sex_erp_&species._tvalue 
        ttests_sex_erp_&species._pval 
        eqvr_sex_erp_&species._val 
        eqvr_sex_erp_&species._pval ;
by erp_plus geneID;
run ;

data &species._results2_erp ;
retain erp_plus geneID;
set &species._results_erp ;
drop _name_ _label_ ;
format pval_ttest_equal best32.30 ;
format pval_ttest_unequal best32.30 ;
run ;

data &species._ttest_results_erp ;
set  &species._results2_erp ;
if pval_ttest_unequal le 0.05 then flag_ttest05_unequal = 1 ;
else  flag_ttest05_unequal = 0 ;
if pval_ttest_equal le 0.05 then flag_ttest05_equal = 1 ;
else  flag_ttest05_equal = 0 ;
run ;
run ;
%mend ;

%ttest (dmel) ;
%ttest (dsim) ;
%ttest (dser) ;


/* merge raw, log uq, oncall flags, ttest for mel, sim and ser */
%macro files (species) ;

proc sort data = &species._erp_wide2 ;
by erp_plus geneID ; ;
proc sort data = &species._uq_erp_wide2 ;
by erp_plus geneID ; ;
proc sort data = sex.onCall_erp_&species._flags ;
by erp_plus geneID ; ;
proc sort data = &species._ttest_results_erp  ;
by erp_plus geneID ; ;
run ;

data dataFile2_&species._erp ;
merge  &species._erp_wide2 &species._uq_erp_wide2 sex.onCall_erp_&species._flags &species._results_erp ;
by erp_plus geneID ; ;
run ;

data dataFile_erp_&species ;
retain erp_plus geneID  ;
set dataFile2_&species._erp  ;
if flag_analyzable = 0 then do ;
    effect_size_equal = .;
    effect_size_unequal = .;
    pval_ttest_equal = . ;
    pval_ttest_unequal = . ;
    probF_folded_f = . ;
    pval_folded_f = .  ;
end ;
drop _name_ _label_ ;
run;
%mend ;

%files (dmel) ;
%files (dsim) ;
%files (dser) ;

/* merge raw, and oncall flags for yak and san */
%macro file2 (species) ;

proc sort data = &species._erp_wide2 ;
by erp_plus geneID ;
proc sort data = sex.onCall_erp_&species._flags ;
by erp_plus geneID ;
run;

data dataFile2_&species._erp ;
merge  &species._erp_wide2  sex.onCall_erp_&species._flags ;
by erp_plus geneID ;
run ;

data dataFile_erp_&species. ;
retain erp_plus geneID ;
set dataFile2_&species._erp  ;
run ;
%mend ;

%file2 (dyak) ;
%file2 (dsan) ;


/* link in dmel6 geneID and symbol */
%macro dmelLink (species, shrt) ; 

data link_&species._2_dmel ;
set sex.fivespecies_&species._w_geneID_02amm ;
where num_&species._geneID =1 and num_dmel6_geneID = 1 ;
keep &species._geneID dmel6_geneID ;
rename &species._geneID = geneID ;
run;

proc sort data= link_&species._2_dmel ;
by geneID ;
proc sort data= dataFile_erp_&shrt. ; 
by geneID ;
run ;

data dataFile_erp_&shrt._ready ;
merge  dataFile_erp_&shrt. (in=in1) link_&species._2_dmel (in=in2) ;
by geneID ;
if in1 ;
run ;

%mend ;

%dmelLink (dsim2, dsim) ;
%dmelLink (dser1, dser) ;
%dmelLink (dyak2, dyak) ;
%dmelLink (dsan1, dsan) ;


data dataFile_erp_dmel_ready ;
set dataFile_erp_dmel ;
run ;

/* make perm  */
%macro perms (species, shrt) ;

data sex.datafile_erp_&species. ;
set dataFile_erp_&shrt._ready ; 
run ;
%mend ;

run ;
%mend ;

%perms (dmel6, dmel) ;
%perms (dsim2, dsim) ;
%perms (dser1, dser) ;
%perms (dyak2, dyak) ;
%perms (dsan1, dsan) ;


/* fix pval formating, add additional bias and sex limited flags and sexClass */

%macro upd_files (species, shrt) ;

/* fix pvalue */
data upd_datafile_erp_&species.;
set  sex.datafile_erp_&species. ;
format pval_ttest_equal best28.24 ;
format pval_ttest_unequal best28.24;
format pval_folded_f best28.24 ;
run;

data check_&species. ;
set  upd_datafile_erp_&species.; 
if find(pval_ttest_equal, "<") ge 1 ;
run;

/* sex_limited_F/M_rc50 */
data upd2_datafile_erp_&species. ;
set  upd_datafile_erp_&species.; 
sum_F = (erpCnts_&shrt._f_rep1 + erpCnts_&shrt._f_rep2 + erpCnts_&shrt._f_rep3 + erpCnts_&shrt._f_rep4 + erpCnts_&shrt._f_rep5 + erpCnts_&shrt._f_rep6 ) ;
sum_M = (erpCnts_&shrt._m_rep1 + erpCnts_&shrt._m_rep2 + erpCnts_&shrt._m_rep3 + erpCnts_&shrt._m_rep4 + erpCnts_&shrt._m_rep5 + erpCnts_&shrt._m_rep6 ) ;
if (sum_F ge 50 and sum_M = 0) then flag_sex_limited_F_rc50 = 1 ;
else flag_sex_limited_F_rc50 = 0;
if (sum_M ge 50 and sum_F = 0) then flag_sex_limited_M_rc50 = 1 ;
else flag_sex_limited_M_rc50 = 0;
drop sum_F sum_M ;
run;

/* sex bias */
data upd3_datafile_erp_&species. ;
set  upd2_datafile_erp_&species.; 
if pval_ttest_equal = . then do ;
    flag_sex_bias_F = . ;
    flag_sex_bias_M = . ;
end ;
else if (pval_ttest_equal le 0.05 and effect_size_equal > 0) then flag_sex_bias_F = 1; 
else if (pval_ttest_equal le 0.05 and effect_size_equal < 0) then flag_sex_bias_M = 1 ;
run ;

/*  sex class */
data upd4_datafile_erp_&species. ;
set  upd3_datafile_erp_&species.; 
length sexClass $12.;
if flag_sex_limited_F_rc50 = 1  then sexClass = "F_limited" ;
else if flag_sex_limited_M_rc50 = 1 then sexClass = "M_limited" ;

else if flag_analyzable = 0 then sexClass = "unanalyzed" ;

else if flag_sex_bias_F = 1 then sexClass = "F_bias" ;
else if flag_sex_bias_M = 1 then sexClass = "M_bias" ;

else if pval_ttest_equal > 0.05 then sexClass = "unbiased";
else sexClass = "no_ttest" ;
run;

title "&species.";
proc freq data = upd4_datafile_erp_&species. ; 
table sexClass ;
run;
title "";
%mend ;

%upd_files (dmel6, dmel) ;
%upd_files (dsim2, dsim) ;
%upd_files (dser1, dser) ;



/* for yak and san:  flag_sex_limited_F/M_rc50, sex bias flags, sexClass*/

%macro upding (species, shrt) ;

/* flag_sex_limited_F/M_rc50 */
data upd_datafile_erp_&species. ;
set  sex.datafile_erp_&species. ; 
if (erpCnts_&shrt._F_rep1 ge 50 and erpCnts_&shrt._F_rep1 = 0) then flag_sex_limited_F_rc50 = 1 ;
else flag_sex_limited_F_rc50 = 0;
if (erpCnts_&shrt._M_rep1 ge 50 and erpCnts_&shrt._M_rep1 = 0) then flag_sex_limited_M_rc50 = 1 ;
else flag_sex_limited_M_rc50 = 0;
run ;

/* sex bias flags */
data upd2_datafile_erp_&species. ;
set  upd_datafile_erp_&species.; 
flag_sex_bias_F = . ; 
flag_sex_bias_M = . ;
if (erpCnts_&shrt._F_rep1 ge 5 or erpCnts_&shrt._M_rep1 ge 5) then flag_5_reads = 1 ; else flag_5_reads = 0 ;
    if erpCnts_&shrt._F_rep1 ge (10*erpCnts_&shrt._M_rep1) and flag_5_reads = 1 then flag_sex_bias_F = 1 ;
    else if erpCnts_&shrt._M_rep1 ge (10*erpCnts_&shrt._F_rep1) and flag_5_reads = 1  then flag_sex_bias_M = 1 ;
drop flag_5_reads ;
run ;

/* sex class */
data upd4_datafile_erp_&species. ;
set upd2_datafile_erp_&species.; 
if flag_sex_limited_F_rc50 = 1 then sexClass = "F_limited" ;
else if flag_sex_limited_M_rc50 = 1 then sexClass = "M_limited" ;
else if flag_analyzable = 0 then sexClass = "unanalyzed" ;
else if flag_sex_bias_F = 1 then sexClass = "F_bias" ;
else if flag_sex_bias_M = 1 then sexClass = "M_bias" ;
else sexClass = "unanalyzed" ;
run;

title "&species.";
proc freq data = upd4_datafile_erp_&species. ; 
table sexClass ;
run;
title "";
%mend ;

%upding (dsan1, dsan) ;
%upding (dyak2, dyak) ;



/* make perm and export */
/* discussed with LMM and decided to drop the following flags */
%macro exp (species, shrt) ;

data sex.datafile_erp_&species._02amm ;
set upd4_datafile_erp_&species.;
drop flag_Fonly_on_readCnt0
flag_Fonly_on_readCnt5
flag_Fonly_on_readCnt10
flag_Fonly_on_readCnt25
flag_Monly_on_readCnt0
flag_Monly_on_readCnt5
flag_Monly_on_readCnt10
flag_Monly_on_readCnt25
flag_both_sex_off_readCnt0
flag_both_sex_off_readCnt5
flag_both_sex_off_readCnt10
flag_both_sex_off_readCnt25
flag_both_sex_on_readCnt0
flag_both_sex_on_readCnt5
flag_both_sex_on_readCnt10
flag_both_sex_on_readCnt25 ;
run ;

proc export data = sex.datafile_erp_&species._02amm
outfile = "!MCLAB/sex_specific_splicing/submission/supplementary/datafiles/datafile_erp_&species..csv"
dbms = csv replace ;
run ;
%mend ;

%exp (dmel6, dmel) ;
%exp (dsim2, dsim) ;
%exp (dser1, dser) ;
%exp (dyak2, dyak) ;
%exp (dsan1, dsan) ;


proc contents data = sex.datafile_erp_dmel6_02amm ; run;

flag_Fonly_on_readCnt0
flag_Fonly_on_readCnt5
flag_Fonly_on_readCnt10
flag_Fonly_on_readCnt25
flag_Monly_on_readCnt0
flag_Monly_on_readCnt5
flag_Monly_on_readCnt10
flag_Monly_on_readCnt25
flag_both_sex_off_readCnt0
flag_both_sex_off_readCnt5
flag_both_sex_off_readCnt10
flag_both_sex_off_readCnt25
flag_both_sex_on_readCnt0
flag_both_sex_on_readCnt5
flag_both_sex_on_readCnt10
flag_both_sex_on_readCnt25

