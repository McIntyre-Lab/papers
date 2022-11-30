

libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";

libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';



/* 
SIM
roll ttest flags up to gene level

input:  dTemp.sim_ttest_frags_with_anno

output:
    work.sim_manual_model_check
    MCLAB/Dros_PB_ChIP/RNAseq/model_output/sim_manual_model_check.csv

*/


/* above input file already contains FBgn - ready to roll up to gene level */

/* only for obs that went through ttest model or counts are going to be off 
        where: flag_multigene = 0
        Where: flag_sim_m_on0_apn=1 or flag_sim_f_on0_apn=1  */



data sim_model_gene2 ;
retain featureID fbgn ;
set dTemp.sim_frag_flags_kitchen_sink ;
where flag_multigene = 0;
if flag_sim_m_on0_apn=1 or flag_sim_f_on0_apn=1 ;
run; /* now have 103921 obs */

proc freq data = sim_model_gene2 noprint ;
tables FBgn / out = cnt_genesagain;
run ;  /* 12598 genes */

proc sort data = sim_model_gene2 ;
by featureType ;
run;

/* create flags for combination of flag_ttest_pval and ratio_expressed */

/* roll flag comb to gene level */

%macro rollup (var, var2, flag1, flag2) ;

proc means data = sim_model_gene2 sum nway noprint ;
class FBgn;
by featureType ;
var &var. ;
output out = sim_combo_&flag1._&var2._&flag2. (drop = _:) sum = num_frags_ttest_&flag1._&var2._&flag2.;
run ;

proc sort data = sim_combo_&flag1._&var2._&flag2.;
by FBgn featureType;
run;

%mend ;

%rollup (flag_ttest_1_ratio2_M , ratio2, 1, M );
%rollup (flag_ttest_1_ratio2_F , ratio2, 1, F );
%rollup (flag_ttest_1_ratio2_U , ratio2, 1, U );

%rollup (flag_ttest_0_ratio2_M , ratio2, 0, M );
%rollup (flag_ttest_0_ratio2_F , ratio2, 0, F );
%rollup (flag_ttest_0_ratio2_U , ratio2, 0, U );

%rollup (flag_ttest_1_trend_M , trend, 1, M );
%rollup (flag_ttest_1_trend_F , trend, 1, F );
%rollup (flag_ttest_1_trend_U , trend, 1, U );

%rollup (flag_ttest_0_trend_M , trend, 0, M );
%rollup (flag_ttest_0_trend_F , trend, 0, F );
%rollup (flag_ttest_0_trend_U , trend, 0, U );

data sim_combo_flags ;
merge sim_combo_: ;
by FBgn featureType ;
run ;

/*for gene level sum the Make flag_manual_model_check  to gene variable = manual_model_check*/
proc sort data = sim_model_gene2 ;
by featureType fbgn ;

proc means data = sim_model_gene2 sum nway noprint ;
class FBgn ;
by featureType ;
var flag_manual_model_check ;
output out = sim_manual_model_check (drop =_:)  sum = manual_model_check;
run;

proc means data = sim_model_gene2 sum  nway noprint ;
class FBgn ;
by featureType ;
var flag_ttest_pval ;
output out = sim_sum_frags_le05 (drop =_:) sum = sum_frags_le05
run;

proc means data = sim_model_gene2 min  nway noprint ;
class FBgn ;
by featureType ;
var ttest_pval ;
output out = sim_ttest_minpval (drop = _:) min = ttest_minpval;
run;

proc freq data = sim_model_gene2 noprint ;
by featureType ;
tables FBgn / out = sim_num_frags_ttest (drop = percent rename = (count = num_frags_ttest));
run;

proc sort data = sim_combo_flags ;
by FBgn featureType ;
proc sort data = sim_num_frags_ttest ;
by FBgn featureType ;
proc sort data = sim_ttest_minpval;
by FBgn featureType ;
proc sort data = sim_sum_frags_le05 ;
by FBgn featureType ;
proc sort data = sim_manual_model_check ;
by FBgn featureType ;
run;

data sim2_manual_model_check ;
merge   sim_num_frags_ttest     sim_manual_model_check  sim_sum_frags_le05  sim_ttest_minpval   sim_combo_flags;
by FBgn featureType ;
run ;  /* 43834 obs */

data dTemp.sim_manual_model_check ;
set  sim2_manual_model_check ;
label num_frags_ttest = "num_frags_ttest";
label ttest_minpval = "ttest_minpval" ;
if manual_model_check = . then num_frags_ttest = .;
run;  


/* take the minimum pvalue from the satherwaite pvalues across fragments keep as min_p*/
/* sum satherwaite flags keep as sum_frags_le05*/
/*need to know how many frags there are num_frags_ttest*/

/*final shoud be  geneid, numfrags,
manual_model_check, sum_frags_le05, satt_minp*/



/* export to share for Adalena 
    note that 'all' fragments in starting file are present even if off */

proc export data = dTemp.sim_manual_model_check
outfile = "!MCLAB/Dros_PB_ChIP/RNAseq/model_output/sim_manual_model_check.csv"
dbms = csv replace ;
run ;






