
libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';



/* 
MEL
roll ttest flags up to gene level

input:  dTemp.mel_ttest_frags_with_anno

output:
    work.mel_manual_model_check
    MCLAB/Dros_PB_ChIP/RNAseq/model_output/mel_manual_model_check.csv
*/



/* above input file already contains FBgn - ready to roll up to gene level */

/* only for obs that went through ttest model or counts are going to be off 
        where: flag_multigene = 0
        Where: flag_mel_m_on0_apn=1 or flag_mel_f_on0_apn=1  */

data mel_model_gene2 ;
retain featureID fbgn ;
set dTemp.mel_frag_flags_kitchen_sink ;
where flag_multigene = 0;
if flag_mel_m_on0_apn=1 or flag_mel_f_on0_apn=1 ;
run; /* now have 115097 obs */

proc freq data = mel_model_gene2 ;
tables FBgn / out = cnt_genesagain;
run ;  /* 13382 genes */


/* roll flag comb to gene level */
proc sort data = mel_model_gene2 ;
by featureType fbgn ;
run;

%macro rollup (var, var2, flag1, flag2) ;

proc means data = mel_model_gene2 sum nway noprint ;
class FBgn;
by featureType ;
var &var. ;
output out = mel_combo_&flag1._&var2._&flag2. (drop = _:) sum = num_frags_ttest_&flag1._&var2._&flag2.;
run ;

proc sort data = mel_combo_&flag1._&var2._&flag2.;
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

data mel_combo_flags ;
merge mel_combo_: ;
by FBgn featureType ;
run ;

/*for gene level sum the Make flag_manual_model_check  to gene variable = manual_model_check*/
proc sort data = mel_model_gene2 ;
by featureType fbgn ;

proc means data = mel_model_gene2 sum nway noprint ;
class FBgn ;
by featureType ;
var flag_manual_model_check ;
output out = mel_manual_model_check (drop =_:)  sum = manual_model_check;
run;

proc means data = mel_model_gene2 sum  nway noprint ;
class FBgn ;
by featureType ;
var flag_ttest_pval ;
output out = mel_sum_frags_le05 (drop =_:) sum = sum_frags_le05
run;

proc means data = mel_model_gene2 min  nway noprint ;
class FBgn ;
by featureType ;
var ttest_pval ;
output out = mel_ttest_minpval (drop = _:) min = ttest_minpval;
run;

proc freq data = mel_model_gene2 noprint ;
by featureType ;
tables FBgn / out = mel_num_frags_ttest (drop = percent rename = (count = num_frags_ttest));
run;

proc sort data = mel_combo_flags ;
by FBgn featureType ;
proc sort data = mel_num_frags_ttest ;
by FBgn featureType ;
proc sort data = mel_ttest_minpval;
by FBgn featureType ;
proc sort data = mel_sum_frags_le05 ;
by FBgn featureType ;
proc sort data = mel_manual_model_check ;
by FBgn featureType ;
run;

data mel2_manual_model_check ;
merge   mel_num_frags_ttest mel_manual_model_check  mel_sum_frags_le05  mel_ttest_minpval   mel_combo_flags;
by FBgn featureType ;
run ;  /* 44742 obs */

data dTemp.mel_manual_model_check ;
set  mel2_manual_model_check ;
label num_frags_ttest = "num_frags_ttest";
label ttest_minpval = "ttest_minpval" ;
if manual_model_check = . then num_frags_ttest = .; /* I don't understand this maybe a shortcut*/
run;  


/* take the minimum pvalue from the satherwaite pvalues across fragments keep as min_p*/
/* sum satherwaite flags keep as sum_frags_le05*/
/*need to know how many frags there are num_frags_ttest*/

/*final shoud be  geneid, numfrags,
manual_model_check, sum_frags_le05, satt_minp combo_flags

note to check that the number_of fragments here when merging to the previous gene level annotations should be the same or smaller than the number of fragments for the number of flags for the gene_ratio2-at the gene level  */


/* export to share for Adalena 
    note that 'all' fragments in starting file are present even if off */

proc export data = dTemp.mel_manual_model_check
outfile = "!MCLAB/Dros_PB_ChIP/RNAseq/model_output/mel_manual_model_check.csv"
dbms = csv replace ;
run ;



