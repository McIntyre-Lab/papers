
/*

all jxnHash data files

*/

libname melsas "!MCLAB/useful_dmel_data/flybase650/sas_data" ;

libname sex "!MCLAB/sex_specific_splicing/sasdata" ;


/* wide format file with raw counts */
%macro raw (species) ;
proc sort data = sex.&species._ujc_cnts_sumtr_stack_w0s ;
by jxnHash ;
run;

proc transpose data = sex.&species._ujc_cnts_sumtr_stack_w0s out = &species._raw_wide prefix = rawCnts_;
by jxnHash ;
id sample ;
var jxnHash_sumReadNum ;
run ;

data &species._raw_wide2 ;
set &species._raw_wide  ;
drop _name_ ;
run;
%mend ;

%raw (dmel) ;
%raw (dsim) ;
%raw (dser) ;
%raw (dyak) ;
%raw (dsan) ;

/* create wide format file with log_uq_sumReadNum for merging in  */
%macro uqs (species, genome) ;

proc sort data= sex.norm_ujc_&species. ;
by jxnHash ;
run;

proc transpose data = sex.norm_ujc_&species. out = &species._uq_ujc_wide prefix = logUQ_;
by jxnHash ;
id sample ;
var log_uq_sumReadNum ;
run;

data &species._uq_ujc_wide2 ;
set &species._uq_ujc_wide;
drop _name_ ;
run ;
%mend ;

%uqs (dmel, dmel6) ;
%uqs (dsim, dsim2) ;
%uqs (dser, dser1) ;


%macro onoffs (species ) ;
/* add on off flags - includes flag analyzable */
proc sort data = sex.onCall_ujc_&species._flags ;
by jxnHash  ;
run ;

%mend ;

%onoffs (dmel) ;
%onoffs (dsim) ;
%onoffs (dser) ;
%onoffs (dyak) ;
%onoffs (dsan) ;

/* ttest results */
%macro ttest (species) ;

proc transpose data = sex.ttests_sex_ujc_&species. out = ttests_sex_ujc_&species._tvalue prefix=effect_size_;
by jxnHash ;
id variances ;
var tValue ;
run ;

proc transpose data = sex.ttests_sex_ujc_&species. out = ttests_sex_ujc_&species._pval prefix=pval_ttest_;
by jxnHash ;
id variances ;
var probt ;
run ;

proc transpose data = sex.eqvr_sex_ujc_&species. out = eqvr_sex_ujc_&species._val prefix=probF_;
by jxnHash ;
id method ;
var Fvalue ;
run ;

proc transpose data = sex.eqvr_sex_ujc_&species. out = eqvr_sex_ujc_&species._pval prefix=pval_;
by jxnHash ;
id method ;
var ProbF ;
run ;

proc sort data= ttests_sex_ujc_&species._tvalue ;
by jxnHash ;
proc sort data= ttests_sex_ujc_&species._pval ;
by jxnHash ;
proc sort data= eqvr_sex_ujc_&species._val ;
by jxnHash ;
proc sort data= eqvr_sex_ujc_&species._pval ;
by jxnHash ;
run ;

data &species._results ;
merge   ttests_sex_ujc_&species._tvalue 
        ttests_sex_ujc_&species._pval 
        eqvr_sex_ujc_&species._val 
        eqvr_sex_ujc_&species._pval ;
by jxnHash;
run ;

data &species._results2 ;
retain jxnHash;
set &species._results ;
drop _name_ _label_ ;
format pval_ttest_equal best32.30 ;
format pval_ttest_unequal best32.30 ;
run ;

data &species._ttest_results ;
set  &species._results2 ;
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


/* merge raw, log uq, oncall flags, ttest for mel, sim and ser [only reps 4, 5, and 6]*/
%macro files (species) ;

proc sort data = &species._raw_wide2 ;
by jxnHash ;
proc sort data = &species._uq_ujc_wide2 ;
by jxnHash ;
proc sort data = sex.onCall_ujc_&species._flags ;
by jxnHash ;
proc sort data = &species._ttest_results  ;
by jxnHash ;
run ;

data dataFile2_&species ;
merge  &species._raw_wide2 sex.onCall_ujc_&species._flags &species._uq_ujc_wide2  &species._results ;
by jxnHash ;
run ;

data dataFile_&species ;
retain jxnHash geneID ;
set dataFile2_&species  ;
if flag_analyzable = 0 then do ;
    effect_size_equal = .;
    effect_size_unequal = .;
    pval_ttest_equal = . ;
    pval_ttest_unequal = . ;
    probF_folded_f = . ;
    pval_folded_f = .  ;
end ;
drop _label_ _name_ ;
run;
%mend ;

%files (dmel) ;
%files (dsim) ;
%files (dser) ;

/* merge raw, and oncall flags for yak and san */
%macro file2 (species) ;

proc sort data = &species._raw_wide2 ;
by jxnHash ;
proc sort data = sex.onCall_ujc_&species._flags ;
by jxnHash ;
run;

data dataFile2_&species ;
merge  &species._raw_wide2  sex.onCall_ujc_&species._flags ;
by jxnHash ;
run ;

data dataFile_&species ;
retain jxnHash geneID ;
set dataFile2_&species  ;
run ;
%mend ;

%file2 (dyak) ;
%file2 (dsan) ;

/* fivespecies genome w geneid 02amm to link jxnHash to dmel jxnHash - sim, ser, yak and san */
%macro links (species, shrt) ;

data &shrt._link ;
set sexsplic.fivespecies_&species._w_geneid_02amm;
rename &species._jxnHash = jxnHash ; 
keep &species._jxnHash dmel6_jxnHash dmel6_geneID num_:;
run ;

proc sort data= &shrt._link nodups ;
by _all_ ;
run;

data almost_&shrt._data ;
merge datafile_&shrt (in=in1) &shrt._link (in=in2) ;
by jxnHash ;
if in1 ;
run ;

%mend ;

%links (dsim2, dsim) ;
%links (dser1, dser) ;
%links (dyak2, dyak) ;
%links (dsan1, dsan) ;

/* make dmel match */
data almost_dmel_data ;
set dataFile_dmel ;
run ;

/* make perm and export */

%macro perming (species, shrt) ;

data sex.datafile_jxnHash_&species. ;
set almost_&shrt._data ;
run ;
%mend ;

%perming (dmel6, dmel) ;
%perming (dsim2, dsim) ;
%perming (dser1, dser) ;
%perming (dyak2, dyak) ;
%perming (dsan1, dsan) ;



/* fix pval formating, add additional bias and sex limited flags and sexClass */

%macro upd_files (species, shrt) ;

/* fix pvalue */
data upd_datafile_jxnHash_&species.;
set  sex.datafile_jxnHash_&species. ;
format pval_ttest_equal best28.24 ;
format pval_ttest_unequal best28.24;
format pval_folded_f best28.24 ;
run;

data check_&species._jxnHash ;
set  upd_datafile_jxnHash_&species.; 
if find(pval_ttest_equal, "<") ge 1 ;
run;

/* sex_limited_F/M_rc50 */
data upd2_datafile_jxnHash_&species. ;
set  upd_datafile_jxnHash_&species.; 
sum_F = (rawCnts_&shrt._f_rep1 + rawCnts_&shrt._f_rep2 + rawCnts_&shrt._f_rep3 + rawCnts_&shrt._f_rep4 + rawCnts_&shrt._f_rep5 + rawCnts_&shrt._f_rep6 ) ;
sum_M = (rawCnts_&shrt._m_rep1 + rawCnts_&shrt._m_rep2 + rawCnts_&shrt._m_rep3 + rawCnts_&shrt._m_rep4 + rawCnts_&shrt._m_rep5 + rawCnts_&shrt._m_rep6 ) ;
if (sum_F ge 50 and sum_M = 0) then flag_sex_limited_F_rc50 = 1 ;
else flag_sex_limited_F_rc50 = 0;
if (sum_M ge 50 and sum_F = 0) then flag_sex_limited_M_rc50 = 1 ;
else flag_sex_limited_M_rc50 = 0;
drop sum_F sum_M ;
run;

/* sex bias */
data upd3_datafile_jxnHash_&species. ;
set  upd2_datafile_jxnHash_&species.; 
if pval_ttest_equal = . then do ;
    flag_sex_bias_F = . ;
    flag_sex_bias_M = . ;
end ;
else if (pval_ttest_equal le 0.05 and effect_size_equal > 0) then flag_sex_bias_F = 1; 
else if (pval_ttest_equal le 0.05 and effect_size_equal < 0) then flag_sex_bias_M = 1 ;
run ;

/*  sex class */
data upd4_datafile_jxnHash_&species. ;
set  upd3_datafile_jxnHash_&species.; 
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
proc freq data = upd4_datafile_jxnHash_&species. ; 
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
data upd_datafile_jxnHash_&species. ;
set  sex.datafile_jxnHash_&species. ; 
where geneID ne "";
if (rawCnts_&shrt._F_rep1 ge 50 and rawCnts_&shrt._F_rep1 = 0) then flag_sex_limited_F_rc50 = 1 ;
else flag_sex_limited_F_rc50 = 0;
if (rawCnts_&shrt._M_rep1 ge 50 and rawCnts_&shrt._M_rep1 = 0) then flag_sex_limited_M_rc50 = 1 ;
else flag_sex_limited_M_rc50 = 0;
run ;

/* sex bias flags */
data upd2_datafile_jxnHash_&species. ;
set  upd_datafile_jxnHash_&species.; 
flag_sex_bias_F = . ; 
flag_sex_bias_M = . ;
if (rawCnts_&shrt._F_rep1 ge 5 or rawCnts_&shrt._M_rep1 ge 5) then flag_5_reads = 1 ; else flag_5_reads = 0 ;
    if rawCnts_&shrt._F_rep1 ge (10*rawCnts_&shrt._M_rep1) and flag_5_reads = 1 then flag_sex_bias_F = 1 ;
    else if rawCnts_&shrt._M_rep1 ge (10*rawCnts_&shrt._F_rep1) and flag_5_reads = 1  then flag_sex_bias_M = 1 ;
drop flag_5_reads ;
run ;

/* sex class */
data upd4_datafile_jxnHash_&species. ;
set upd2_datafile_jxnHash_&species.; 
if flag_sex_limited_F_rc50 = 1 then sexClass = "F_limited" ;
else if flag_sex_limited_M_rc50 = 1 then sexClass = "M_limited" ;
else if flag_analyzable = 0 then sexClass = "unanalyzed" ;
else if flag_sex_bias_F = 1 then sexClass = "F_bias" ;
else if flag_sex_bias_M = 1 then sexClass = "M_bias" ;
else sexClass = "unanalyzed" ;
run;

title "&species.";
proc freq data = upd4_datafile_jxnHash_&species. ; 
table sexClass ;
run;
title "";
%mend ;

%upding (dsan1, dsan) ;
%upding (dyak2, dyak) ;



/* drop unneeded flags (discussed which with LMM), make perm and export */
%macro exp (species, shrt) ;

data sex.datafile_jxnHash_&species._02amm ;
set upd4_datafile_jxnHash_&species.;
drop 
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
flag_both_sex_on_readCnt25 ;
run ;

proc export data = sex.datafile_jxnHash_&species._02amm
outfile = "!MCLAB/sex_specific_splicing/Tables/datafile_jxnHash_&species..csv"
dbms = csv replace ;
run ;
%mend ;

%exp (dmel6, dmel) ;
%exp (dsim2, dsim) ;
%exp (dser1, dser) ;
%exp (dyak2, dyak) ;
%exp (dsan1, dsan) ;


proc contents data = sex.datafile_jxnHash_dmel6_02amm ; run;

/* drop following flags */



