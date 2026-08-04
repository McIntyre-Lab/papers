
libname lmm_rmg "!MCLAB/lmm_rmg_dros_data/sasdata";
libname lmm_rlr "!MCLAB/lmm_rlr_head_data/sasdata";
libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata";
libname sex "!MCLAB/sex_specific_splicing/sasdata";

/*


flag on/off in ERP count files
using ERP_plus var 

* if ge 50% of reps then is expressed

using readNums of 0, 5, 10 and 25:  
    if off in F and off in M then   flag_both_sex_off_readNum&num
    if on in F and on on M then     flag_both_sex_on_readNum&num
    if off in F and on in M then    flag_Monly_on_readNum&num
    if on in F and off in M then    flag_Fonly_on_readNum&num


for norm:
   drop if flag_both_sex_off_readNum10 = 1 

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

data design_rr ;
set lmm_rlr.design_rlr_yak_san ;
new = compress(sample||"_rep1") ;
drop sample ;
rename new = sample ;
rep = "rep1";
keep new sex rep species ;
run;

proc sort data = design_rr nodups ;
by _all_ ;
run;


%macro import (species, ref, design) ;

proc sort data= &design ;
by sample ;
proc sort data = sex.&species._erp_cnts_sumtr_stack_w0s ;
by sample ;
run ;

data &species._erp_cnts ;
merge &design (in=in1) sex.&species._erp_cnts_sumtr_stack_w0s (in=in2) ; 
by sample ;
if in2 ;
run;

data &species._erp_cnts1;
set &species._erp_cnts ;
if ERP_sumReadNum > 0 then readCnt_on0 = 1; else readCnt_on0 = 0;
if ERP_sumReadNum >= 5 then readCnt_on5 = 1; else readCnt_on5 = 0;
if ERP_sumReadNum >= 10 then readCnt_on10 = 1; else readCnt_on10 = 0;
if ERP_sumReadNum >= 25 then readCnt_on25 = 1; else readCnt_on25 = 0;
run;

data &species._erp_cnts2 ;
set &species._erp_cnts1 ;
if (species = "dmel" or species = "dsim" or species = "dser") and (rep = "rep4" or rep = "rep5" or rep = "rep6") then do ;
output ;
end ;
else output ;
run ; 

%mend ;

%import (dmel, dmel6, design_r) ; 
%import (dsim, dsim2, design_r) ; 
%import (dser, dser1, design_k) ; 
%import (dyak, dyak2, design_rr) ; 
%import (dsan, dsan1, design_rr) ; 


%macro which (species, ref) ;

%macro oncalls (num) ;

/* Going to use readNum 0, 5, 10 and 25 for flagging on/off */
proc sort data= &species._erp_cnts2 ;
by sex geneID ERP_plus ;
run;

/* gene is on if expressed in ge 50% of the reps */
proc means data = &species._erp_cnts2  noprint ;
    by sex geneID ERP_plus  ;
    var readCnt_on&num;
    output out = &species._sex_on&num mean=sex_percent_on ;
    run ;

proc sort data = &species._sex_on&num ;
    by geneID ERP_plus ;
    run;

proc transpose data = &species._sex_on&num out = &species._sex_on&num._sbys ;
    by geneID ERP_plus ;
    id sex ;
    var sex_percent_on ;
    run;

               * if 50% of reps then is expressed ;
               * using sex to determine ;


data sex.onCall_ERP_&species._on&num.;  
set &species._sex_on&num._sbys ;   
    if  f ge 0.5 then flag_F_on = 1 ; else flag_F_on = 0 ;
    if  m ge 0.5 then flag_M_on = 1 ; else flag_M_on = 0 ;

    /* if off in F and off in M then flag all off */
    if flag_F_on = 0 and flag_M_on = 0 then flag_both_sex_off_readCnt&num = 1 ;
    else flag_both_sex_off_readCnt&num = 0;

    /* if on in F and on in M then flag all on */
    if flag_F_on = 1 and flag_M_on = 1 then flag_both_sex_on_readCnt&num = 1 ;
    else flag_both_sex_on_readCnt&num = 0 ;

    /* if off in F and on in M then flag M only on */
    if (flag_F_on = 0 and flag_M_on = 1) then flag_Monly_on_readCnt&num  = 1;
    else flag_Monly_on_readCnt&num = 0 ;

    /* if on in F and off in M then flag F only on */
    if flag_F_on = 1 and flag_M_on = 0 then flag_Fonly_on_readCnt&num = 1 ;
    else flag_Fonly_on_readCnt&num = 0 ;

    /* sex limited flags */
    if f > 0 then flag_F_detect_rc&num = 1 ; else flag_F_detect_rc&num = 0 ;
    if m > 0 then flag_M_detect_rc&num = 1 ; else flag_M_detect_rc&num = 0 ;

    if flag_M_detect_rc&num  = 1 and flag_F_detect_rc&num  = 1 then flag_sex_limited_rc&num  = 0 ;
    else if flag_M_detect_rc&num = 0 and flag_F_detect_rc&num  = 0 then flag_sex_limited_rc&num = . ;  /* if F and M both 0 then missing */
    else flag_sex_limited_rc&num = 1 ;
 
    keep geneID ERP_plus flag_: f m ;
    run ;

proc sort data = sex.oncall_erp_&species._on&num. ;
by geneID ERP_plus ;
run ;

%mend ;

%onCalls (0) ;
%onCalls (5) ;
%onCalls (10) ;
%onCalls (25) ;
%mend ;
%which (dmel, dmel6) ; 
%which (dsim, dsim2) ;
%which (dser, dser1) ;

%which (dyak, dyak2) ;
%which (dsan, dsan1) ;

/* create dataset with: 
flag_sex_limited_rc0
flag_both_sex_off_readCnt0, 5, 10, 25 
flag_both_sex_on_readCnt&num 
flag_Monly_on_readCnt&num
flag_Fonly_on_readCnt&num
*/


%macro flags (species) ;

data onCall_ERP_&species._flags2;
merge sex.oncall_erp_&species._on: ;
by geneID ERP_plus ;
run;

data onCall_ERP_&species._flags;
set onCall_ERP_&species._flags2;
keep geneID ERP_plus 
    flag_sex_limited_rc0 flag_both_sex_off_readCnt: flag_both_sex_on_readCnt: flag_Monly_on_readCnt: flag_Fonly_on_readCnt: ;
run ;

data sex.onCall_ERP_&species._flags; 
retain geneID ERP_plus flag_sex_limited_rc0 flag_analyzable flag_both_sex_off_readCnt0 flag_both_sex_off_readCnt5 flag_both_sex_off_readCnt10 flag_both_sex_off_readCnt25 
flag_both_sex_on_readCnt0 flag_both_sex_on_readCnt5 flag_both_sex_on_readCnt10 flag_both_sex_on_readCnt25
flag_Monly_on_readCnt0 flag_Monly_on_readCnt5 flag_Monly_on_readCnt10 flag_Monly_on_readCnt25
flag_Fonly_on_readCnt0 flag_Fonly_on_readCnt5 flag_Fonly_on_readCnt10 flag_Fonly_on_readCnt25 ;
set onCall_ERP_&species._flags; 

if flag_sex_limited_rc0 = 0 and flag_both_sex_off_readCnt5 = 0 then flag_analyzable = 1 ;
else flag_analyzable = 0 ;


run ;
%mend ;

%flags (dmel) ;
%flags (dsim) ;
%flags (dser) ;

%flags (dyak) ;
%flags (dsan) ;

ods pdf file = "!MCLAB/sex_specific_splicing/expression_output/freqs_onCalls_ERP_counts.pdf" ;

%macro which (species) ;
%macro freqs (num) ;

title "&species. MF on calls for ERP_plus,  readCnt > &num";
proc freq data = sex.onCall_ERP_&species._on&num. ;
tables flag_:  ;
run;

proc freq data =sex.onCall_ERP_&species._on&num.;
tables flag_F_on * flag_M_on  ;
run;

title "";
%mend ;


%freqs (0) ;
%freqs (5) ;
%freqs (10) ;
%freqs (25) ;
%mend ;
%which (dmel) ; 
%which (dsim) ;
%which (dser) ;
%which (dyak) ;
%which (dsan) ;

ods pdf close ;


proc freq data = sex.onCall_erp_dmel_flags ;
tables flag_sex_limited_rc0 * flag_both_sex_off_readCnt5  ;
run;
proc freq data = sex.onCall_erp_dmel_flags ;
tables flag_sex_limited_rc0  flag_both_sex_off_readCnt25  ;
run;

/* 530 ERPs that have at least 5 reads in 1 sex and 0 reads in opposite sex 
   of these, there are 14 that have at least 25 reads in 1 sex and 0 reads in opposite sex  (modify above to flag_both_sex_off_readCnt25)

   54873 ERPs that are analyzable (not sex-limited and have at least 5 reads in ge 50% of the reps in at least 1 sex
*/


ods pdf file = "!MCLAB/sex_specific_splicing/expression_output/freqs_onCalls_ERP_sexLimited_flagAnalyze.pdf" ;

%macro which (species) ;

title "&species. MF sex limited counts";
proc freq data = sex.onCall_erp_&species._flags ;
tables flag_sex_limited_rc0 * flag_both_sex_off_readCnt5  ;
run;

title "&species. MF flag_analyze";
proc freq data = sex.onCall_erp_&species._flags ;
tables flag_analyzable  ;
run;

title "";

%mend ;

%which (dmel) ;
%which (dsim) ;
%which (dser) ;
%which (dyak) ;
%which (dsan) ;

ods pdf close ;


