
/* libname dros 'C:\a1stuff\dros';*/

libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';

libname fb "!MCLAB/useful_dmel_data/flybase617/sas_data"; 


/*

input dtemp.mel_ready_4_anova created in:    import_mel_rna_4_anova_02amm.sas

*/


proc sort data=dTemp.mel_ready_4_anova;
by featureid;
run ;

data dtemp.mel;
set dTemp.mel_ready_4_anova;
obs=_n_;
log_apn_uq_ff=log(apn_uq_ff);
if log_apn_uq_ff=. then log_apn_uq_ff=0;  /*think more here*/
run;


proc glimmix data = dtemp.mel plots = pearsonpanel;
class sex featureid;
model log_apn_uq_ff=sex ;
run;

/*this overall model does not look bad for residuals*/

ods pdf file="!MCLAB/Dros_PB_ChIP/RNAseq/model_output/mel_mixed_class_sex_flag_resid.pdf"; 

ods exclude all;
proc mixed data=dtemp.mel;
by featureid;
class sex;
model log_apn_uq_ff=sex / residual outp=dTemp.mel_resid;
repeated /group=sex;
lsmeans sex/diff;
ods output 
	lsmeans = mel_lsmean
	diffs=mel_diffs;
run;
quit;

ods exclude none;
ods results;

proc univariate data=dTemp.mel_resid normal noprint;
by featureid;
var studentresid;
output out=mel_student probn=probn;
run;

data dtemp.mel_student_flag;
set mel_student;
if probn le 0 then flag_resid=.;
else if probn le 0.05 then flag_resid=1;
else flag_resid=0;
run;

proc freq data=dtemp.mel_student_flag;
tables flag_resid;
run;
ods pdf ;

*16.3% is high...;
ods pdf close ; 


ods exclude all;

proc ttest data= dtemp.mel ; 
by featureid;
class sex;
var log_apn_uq_ff;
ods output  equality=dTemp.mel_equality ttests=dTemp.mel_ttests ;
run;

ods exclude none;
ods results;




