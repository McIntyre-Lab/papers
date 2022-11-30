
libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';


/*libname fb "!MCLAB/useful_dmel_data/flybase617/sas_data";  */

/*
MEL
create fragment flags from ttest model 4 DE sex

input 
    ttest results from mel_frags_model_4_de_sex.sas
        dTemp.mel_equality
        dtemp.mel_ttests
        dtemp.mel_student_flag

output
    dTemp.mel_ttest_frags_with_anno
    MCLAB/Dros_PB_ChIP/RNAseq/model_output/mel_ttest_frags_with_anno.csv


*/


/*save the equality test pvalue and flag this is levine's test
for equality of varaince. Keep featureid pvalue and make a flag if it is le 0.05 called flag_levine*/

data levine ;
set dTemp.mel_equality ;
rename probF = levine_pval ;
label probF = "levine_pval" ;
if probF le 0.05 then flag_levine_pval = 1;
else flag_levine_pval = 0 ;
keep featureID probF flag_levine_pval;
run;

/*split ttests into two files one where method= 'pooled" and one "satherwaite" 
keep featureID and pvalues  merge by feature id
make two flags one for pooled one for satherwaite at 0.05... 
is there a differnce in decision?*/

data pool ;
set dtemp.mel_ttests ;
where method = "Pooled" ;
rename probt = ttest_pval ;
label probt = "ttest_pval" ;
if probt le 0.05 then flag_ttest_pval = 1;
else flag_ttest_pval = 0;
keep probt featureid flag_ttest_pval ;
run ;

data satterthwaite ;
set dtemp.mel_ttests ;
where method = "Satterthwaite" ;
rename probt = satterthwaite_pval ;
label probt = "satterthwaite_pval" ;
if probt le 0.05 then flag_satterthwaite_pval =1 ;
else flag_satterthwaite_pval = 0;
keep probt featureid flag_satterthwaite_pval ;
run;

proc sort data = pool ;
by featureid ;
proc sort data = satterthwaite ;
by featureid ;
run;

data both ;
merge pool (in=in1) satterthwaite (in=in2) ;
by featureid ;
if in1 and in2 then output both ;
run;

proc freq data = both ;
tables flag_ttest_pval * flag_satterthwaite_pval ;
run;
/*
flag_ttest_pval
          flag_satterthwaite_pval

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 117282 |      0 | 117282
         |  80.51 |   0.00 |  80.51
         | 100.00 |   0.00 |
         |  99.60 |   0.00 |
---------+--------+--------+
       1 |    467 |  27924 |  28391
         |   0.32 |  19.17 |  19.49
         |   1.64 |  98.36 |
         |   0.40 | 100.00 |
---------+--------+--------+
Total      117749    27924   145673
            80.83    19.17   100.00
*/



/*assuming "no" differnce choose pooled (default)... 
    not a "huge" difference (148 total outof 52790 (148 sig pool only and 0 sig satter only)--> amm choosing pooled */


/*
stick to pooled (ttest)

 flag_manual_model_check = 1
	if levine pval le .05 OR residual pval le 0.05   */

data resids ;
set dtemp.mel_student_flag ;
label flag_resid = "flag_resid_pval";
label probn = "resid_pval" ;
rename flag_resid = flag_resid_pval ;
rename probn = resid_pval ;
run ;

proc sort data = resids ;
by featureid ;
proc sort data = levine ;
by featureid ;
run;

data check ;
merge resids (in=in1) levine (in=in2) ;
by featureid ;
run;

data model_check ;
set check ;
if flag_levine_pval = 1 or flag_resid_pval = 1 then flag_manual_model_check = 1 ;
else flag_manual_model_check =0 ;
run ;

proc freq data = model_check ;
tables flag_manual_model_check ;
run;  /* 44696 obs where flag_manual_model_check = 1 */

proc freq data = model_check ;
tables flag_levine_pval * flag_resid_pval ;
run ;  /* 44 obs missing resid */

/* merge in ttest pval and flag */
proc sort data = model_check ;
by featureID ;
proc sort data = pool ;
by featureID ;
run;

data dTemp.model_pool_mel oops;
merge model_check (in=in1) pool (in=in2) ;
by featureID ;
if in1 and in2 then output dTemp.model_pool_mel ;
else output oops ; /* 0 in oops */
run;

/* want to output frag ttest results with anno flags */

data mel_ids ;
set chiprna.mel_chip_rna_frag_flags_anno ;
run;

proc sort data = mel_ids ;
by featureID ;
proc sort data = dTemp.model_pool_mel ;
by featureID ;
run;

data mel_frags_w_anno  ;
merge dTemp.model_pool_mel (in=in1)  mel_ids (in=in2) ;
by featureID ;
run ;

data mel_ttest_frags_with_anno;
retain featureID FBgn flag_multigene flag_mel_m_on0_apn flag_mel_f_on0_apn ;
set mel_frags_w_anno ;
run;  /* 198811 obs */

/* create combination of flag_ttest_pval and ratio_expressed and ratio_trend*/
proc contents data = mel_ttest_frags_with_anno ; run;

data mel_frag_flags_kitchen_sink;
set mel_ttest_frags_with_anno ;


if flag_ttest_pval = 1 and ratio_expressed ="male" then flag_ttest_1_ratio2_M = 1 ;
	else flag_ttest_1_ratio2_M = 0;
if flag_ttest_pval = 1 and ratio_expressed ="fem" then flag_ttest_1_ratio2_F = 1 ;
	else flag_ttest_1_ratio2_F = 0;
if flag_ttest_pval = 1 and ratio_expressed ="unb" then flag_ttest_1_ratio2_U = 1 ;
	else flag_ttest_1_ratio2_U= 0;

if flag_ttest_pval = 0 and ratio_expressed ="male" then flag_ttest_0_ratio2_M = 1 ;
	else flag_ttest_0_ratio2_M = 0;
if flag_ttest_pval = 0 and ratio_expressed ="fem" then flag_ttest_0_ratio2_F = 1 ;
	else flag_ttest_0_ratio2_F = 0;
if flag_ttest_pval = 0 and ratio_expressed ="unb" then flag_ttest_0_ratio2_U = 1 ;
	else flag_ttest_0_ratio2_U= 0;


if flag_ttest_pval = 1 and ratio_trend ="male" then flag_ttest_1_trend_M = 1 ;
	else flag_ttest_1_trend_M = 0;
if flag_ttest_pval = 1 and ratio_trend ="fem" then flag_ttest_1_trend_F = 1 ;
	else flag_ttest_1_trend_F = 0;
if flag_ttest_pval = 1 and ratio_trend ="unb" then flag_ttest_1_trend_U = 1 ;
	else flag_ttest_1_trend_U= 0;

if flag_ttest_pval = 0 and ratio_trend ="male" then flag_ttest_0_trend_M = 1 ;
	else flag_ttest_0_trend_M = 0;
if flag_ttest_pval = 0 and ratio_trend ="fem" then flag_ttest_0_trend_F = 1 ;
	else flag_ttest_0_trend_F = 0;
if flag_ttest_pval = 0 and ratio_trend ="unb" then flag_ttest_0_trend_U = 1 ;
	else flag_ttest_0_trend_U= 0;

if flag_ttest_pval = . then do 
					flag_ttest_1_ratio2_M = . ;
					flag_ttest_1_trend_M = .;
					flag_ttest_0_ratio2_M = . ;
					flag_ttest_0_trend_M = .;

					end;
if flag_ttest_pval = . then do
					flag_ttest_1_ratio2_F = . ;
					flag_ttest_1_trend_F = . ;
					flag_ttest_0_ratio2_F = . ;
					flag_ttest_0_trend_F = . ;
					end;
if flag_ttest_pval = . then do
					flag_ttest_1_ratio2_U = . ;
					flag_ttest_1_trend_U = . ;
					flag_ttest_0_ratio2_U = . ;
					flag_ttest_0_trend_U = . ;
					end;

run ;

/*check these new flags ratio2 male should be a subset of trend  --> yes, this is true*/
%macro checks (var1) ;
%macro checks2 (var2) ;

proc freq data = mel_frag_flags_kitchen_sink ;
tables flag_ttest_&var1._ratio2_&var2. * flag_ttest_&var1._trend_&var2. / out = check_&var1._&var2. ;
run ;

%mend ;
%checks2 (M);
%checks2 (F) ;
%mend ;
%checks (0) ;
%checks (1) ;

data dTemp.mel_frag_flags_kitchen_sink ;
set mel_frag_flags_kitchen_sink ;
run ;


/* export  */
proc export data = dTemp.mel_frag_flags_kitchen_sink
outfile = "!MCLAB/Dros_PB_ChIP/RNAseq/model_output/mel_frag_flags_kitchen_sink.csv"
dbms = csv replace ;
run ;

