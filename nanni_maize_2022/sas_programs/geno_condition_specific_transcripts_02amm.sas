
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";


/*

are there genes / transcripts detected in only 1 genotype x condition?


yes - but the overall mean for these is~ 100X less than the overall mean for transcripts detected in all 10 genotype x condition


for transcript NOT going into tappas: mean range:  0.067 to 1.169 
for transcript GOING into tappas: mean range:  65.888 to 65.996 

*/



data try ;
set pacbio.rsem_exp_matrix_tpm0_w_flags ;
run ;

proc contents data = try ; run;

proc freq data = try noprint;
tables flag_B73_amb_on0 *flag_B73_ele_on0 * 
flag_C123_amb_on0 *flag_C123_ele_on0 * 
flag_Hp301_amb_on0 *flag_Hp301_ele_on0 * 
flag_Mo17_amb_on0 *flag_Mo17_ele_on0 * 
flag_NC338_amb_on0 *flag_NC338_ele_on0 / out = cnts ;
run ;

data sums ;
set cnts ;
sum = sum(flag_B73_amb_on0 , flag_B73_ele_on0, flag_C123_amb_on0 , flag_C123_ele_on0 ,flag_Hp301_amb_on0, flag_Hp301_ele_on0, flag_Mo17_amb_on0, flag_Mo17_ele_on0, 
flag_NC338_amb_on0, flag_NC338_ele_on0) ;
run;


proc sort data = sums ;
by  sum ;
run;


%macro flags (geno, trt, val ) ;

data uniq_&geno._&trt._&val. ;
set pacbio.rsem_exp_matrix_tpm0_w_flags ;
sum = sum(flag_B73_amb_on0 , flag_B73_ele_on0, flag_C123_amb_on0 , flag_C123_ele_on0 ,flag_Hp301_amb_on0, flag_Hp301_ele_on0, flag_Mo17_amb_on0, flag_Mo17_ele_on0, 
flag_NC338_amb_on0, flag_NC338_ele_on0) ;
if sum = &val. and flag_&geno._&trt._on0 =1 then output uniq_&geno._&trt._&val. ;
run;

proc sort data = uniq_&geno._&trt._&val. ;
by transcriptID ;

proc transpose data = uniq_&geno._&trt._&val. out = flip_&geno._&trt._&val. ;
by transcriptID ;
var &geno._: ;
run ;

data flip2_&val._&geno._&trt. ;
set flip_&geno._&trt._&val. ;
rename _name_ = sample ;
rename col1 = TPM ;
run ;

proc sort data = flip2_&val._&geno._&trt. ;
by sample ;

title "sum = &val., trt = &trt. ";
proc means data = flip2_&val._&geno._&trt. ;
var tpm ;
by sample ;
output out = uniq2_&val._means mean = mean ;
run;
%mend ;

%flags (B73, ele, 1) ;
%flags (C123, ele, 1) ;
%flags (Hp301, ele, 1) ;
%flags (Mo17, ele, 1) ;
%flags (NC338, ele, 1) ;

%flags (B73, amb, 1) ;
%flags (C123, amb, 1) ;
%flags (Hp301, amb, 1) ;
%flags (Mo17, amb, 1) ;
%flags (NC338, amb, 1) ;

%flags (B73, ele, 10) ;
%flags (C123, ele, 10) ;
%flags (Hp301, ele, 10) ;
%flags (Mo17, ele, 10) ;
%flags (NC338, ele, 10) ;

%flags (B73, amb, 10) ;
%flags (C123, amb, 10) ;
%flags (Hp301, amb, 10) ;
%flags (Mo17, amb, 10) ;
%flags (NC338, amb, 10) ;

data all_1 ;
format sample $16. ;
set uniq2_1_: ;
run ;

proc sort data = all_1 ;
by tpm ;
run ;


data all_10 ;
format sample $16. ;
set uniq2_10_: ;
run ;

proc sort data = all_10 ;
by tpm ;
run ;


/* expression levels where (1) didn't go to tappas and (2) went into tappas */
data try ;
set pacbio.rsem_exp_matrix_tpm0_w_flags ;
run ;

proc contents data = try ; run;

data lows ;
set try ;
where flag_into_tappas = 0  ;
run ; 

data highs ;
set try ;
where flag_into_tappas = 1  ;
run ; 


%macro flips (geno, datain ) ;

proc transpose data = &datain. out = flip_&datain._&geno. ;
by transcriptID ;
var &geno._: ;
run ;

data flip2_&datain._&geno. ;
set flip_&datain._&geno. ;
rename _name_ = sample ;
rename col1 = TPM ;
run ;

proc sort data = flip2_&datain._&geno.;
by sample ;


title " geno: &geno. ones not into tappas ";
proc means data = flip2_&datain._&geno.;
var tpm ;
by sample ;
output out = means_&datain._&geno. mean = mean;
run;

data means2_&datain._&geno ;
format genotype $6.;
set means_&datain._&geno ;
genotype = "&geno.";
if find(sample,"Amb") ge 1 then trt = "Amb";
else trt = "Ele" ;
drop _type_ _freq_ ;
run;

%mend ;

%flips (B73, lows) ;
%flips (C123, lows) ;
%flips (Hp301, lows) ;
%flips (Mo17, lows) ;
%flips (NC338, lows) ;

%flips (B73, highs) ;
%flips (C123, highs) ;
%flips (Hp301, highs) ;
%flips (Mo17, highs) ;
%flips (NC338, highs) ;


%macro all (datain) ;

data all_&datain. ;
format sample $16.;
set means2_&datain._: ;
run ;

proc sort data = all_&datain. ;
by genotype trt ;
run;

proc means data = all_&datain. ;
var mean ;
by genotype trt ;
output out = geno_by_trt_&datain. min = min max = max;
run ;

proc sort data = mean_geno_by_trt_&datain. ;
by mean ;
run;
%mend ;

%all (lows ) ;
%all (highs ) ;



proc sort data = all_lows ;
by mean ;
run;   /* mean range:  0.067 to 1.169 */

proc sort data = all_highs ;
by mean ;
run;   /* mean range:  65.888 to 65.996 */





