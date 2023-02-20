
libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";


libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';


/*libname fbSim "!MCLAB/useful_dsim_data/flybase202/sas_data"; */

/*
SIM
create fragment flags from ttest model 4 DE sex

input 
    ttest results from sim_frags_model_4_de_sex.sas
        dTemp.sim_equality
        dtemp.sim_ttests
        dtemp.sim_student_flag

output
    dTemp.sim_ttest_frags_with_anno
    MCLAB/Dros_PB_ChIP/RNAseq/model_output/sim_ttest_frags_with_anno.csv

*/





/*save the equality test pvalue and flag this is levine's test
for equality of varaince. Keep featureid pvalue and make a flag if it is le 0.05 called flag_levine*/

data sim_levine ;
set dTemp.sim_equality ;
if probF le 0.05 then flag_levine_pval = 1;
else flag_levine_pval = 0 ;
keep featureID probF flag_levine_pval;
rename probF = levine_pval ;
label probF = "levine_pval" ;
run;

/* pooled is default */

/*split ttests into two files one where method= 'pooled" and one "satherwaite" 
keep featureID and pvalues  merge by feature id
make two flags one for pooled one for satherwaite at 0.05... 
is there a differnce in decision?*/

data sim_pool ;
set dtemp.sim_ttests ;
where method = "Pooled" ;
rename probt = ttest_pval ;
label probt = "ttest_pval" ;
if probt le 0.05 then flag_ttest_pval = 1;
else flag_ttest_pval = 0;
keep probt featureid flag_ttest_pval ;
run ;

data sim_satterthwaite ;
set dtemp.sim_ttests ;
where method = "Satterthwaite" ;
rename probt = satterthwaite_pval ;
label probt = "satterthwaite_pval" ;
if probt le 0.05 then flag_satterthwaite_pval =1 ;
else flag_satterthwaite_pval = 0;
keep probt featureid flag_satterthwaite_pval ;
run;

proc sort data = sim_pool ;
by featureid ;
proc sort data = sim_satterthwaite ;
by featureid ;
run;

data both ;
merge sim_pool (in=in1) sim_satterthwaite (in=in2) ;
by featureid ;
if in1 and in2 then output both ;
run;


proc freq data = both ;
tables flag_ttest_pval * flag_satterthwaite_pval ;
run;
/*flag_ttest_pval
          flag_satterthwaite_pval

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 118593 |    158 | 118751
         |  87.55 |   0.12 |  87.67
         |  99.87 |   0.13 |
         |  99.48 |   0.97 |
---------+--------+--------+
       1 |    624 |  16077 |  16701
         |   0.46 |  11.87 |  12.33
         |   3.74 |  96.26 |
         |   0.52 |  99.03 |
---------+--------+--------+
Total      119217    16235   135452
            88.01    11.99   100.00
  */


/*
stick to pooled, not satterthwaite

 flag_manual_model_check = 1
	if levine pval le .05 OR resid pval le 0.05 resid  */


data sim_resids ;
set dtemp.sim_student_flag ;
label flag_resid = "flag_resid_pval";
label probn = "resid_pval" ;
rename flag_resid = flag_resid_pval ;
rename probn = resid_pval ;
run ;

proc sort data = sim_resids ;
by featureid ;
proc sort data = sim_levine  ;
by featureid ;
run;

data check ;
merge sim_resids (in=in1) sim_levine (in=in2) ;
by featureid ;
run;

data model_check ;
set check ;
if flag_levine_pval = 1 or flag_resid_pval = 1 then flag_manual_model_check = 1 ;
else flag_manual_model_check =0 ;
run ;

proc freq data = model_check ;
tables flag_manual_model_check ;
run;  /* 41251 obs where flag_manual_model_check = 1 */

proc freq data = model_check ;
tables flag_levine_pval * flag_resid_pval ;
run ;  /* 30 obs missing flag_resid_pval */


/* merge in ttest pval and flag */
proc sort data = model_check ;
by featureID ;
proc sort data = sim_pool ;
by featureID ;
run;

data dTemp.model_pool_sim oops;
merge model_check (in=in1) sim_pool (in=in2) ;
by featureID ;
if in1 and in2 then output dTemp.model_pool_sim ;
else output oops ; /* 0 in oops */
run;


/* want to output frag ttest results with anno flags */

data sim_ids ;
set chiprna.sim_chip_rna_frag_flags_anno ;
run;

proc sort data = sim_ids ;
by featureID ;
proc sort data = dTemp.model_pool_sim ;
by featureID ;
run;

data sim_frags_w_anno  ;
merge dTemp.model_pool_sim (in=in1)  sim_ids (in=in2) ;
by featureID ;
run ;

data sim_ttest_frags_with_anno;
retain featureID FBgn flag_multigene flag_sim_m_on0_apn flag_sim_f_on0_apn ;
set sim_frags_w_anno ;
run;  /* 198811 obs */

proc contents data = sim_ttest_frags_with_anno; run;

/* create combination of flag_ttest_pval and ratio_expressed and ratio_trend*/


data sim_frag_flags_kitchen_sink;
set sim_ttest_frags_with_anno ;


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

proc freq data = sim_frag_flags_kitchen_sink ;
tables flag_ttest_&var1._ratio2_&var2. * flag_ttest_&var1._trend_&var2. / out = check_&var1._&var2. ;
run ;

%mend ;
%checks2 (M);
%checks2 (F) ;
%mend ;
%checks (0) ;
%checks (1) ;

data dTemp.sim_frag_flags_kitchen_sink ;
set sim_frag_flags_kitchen_sink ;
run ;


/* export  */
proc export data = dTemp.sim_frag_flags_kitchen_sink
outfile = "!MCLAB/Dros_PB_ChIP/RNAseq/model_output/sim_frag_flags_kitchen_sink.csv"
dbms = csv replace ;
run ;

