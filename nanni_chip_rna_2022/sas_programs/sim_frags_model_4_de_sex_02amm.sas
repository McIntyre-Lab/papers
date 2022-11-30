
libname dros 'C:\a1stuff\dros';

libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';

libname fbSim "!MCLAB/useful_dsim_data/flybase202/sas_data"; 


/*

input work.sim_ready_4_anova created in:    import_sim_rna_4_anova_02amm.sas

*/


proc sort data=dTemp.sim_ready_4_anova;
by featureid;
run ;

data dtemp.sim;
set dTemp.sim_ready_4_anova;
obs=_n_;
log_apn_uq_ff=log(apn_uq_ff);
if log_apn_uq_ff=. then log_apn_uq_ff=0;  /*think more here  */
run;


proc glimmix data = dtemp.sim plots = pearsonpanel;
class sex featureid;
model log_apn_uq_ff=sex ;
run;

/*this overall model does not look bad for residuals*/

ods pdf file="!MCLAB/Dros_PB_ChIP/RNAseq/model_output/sim_frag_mixed_class_sex_flag_resid.pdf"; 

ods exclude all;
proc mixed data=dtemp.sim;
by featureid;
class sex;
model log_apn_uq_ff=sex / residual outp=dTemp.sim_resid;
repeated /group=sex;
lsmeans sex/diff;
ods output 
	lsmeans = sim_lsmean
	diffs=sim_diffs;
run;
quit;

ods exclude none;
ods results;

proc univariate data=dTemp.sim_resid normal noprint;
by featureid;
var studentresid;
output out=sim_student probn=probn;
run;

data dtemp.sim_student_flag;
set sim_student;
if probn le 0 then flag_resid=.;
else if probn le 0.05 then flag_resid=1;
else flag_resid=0;
run;

proc freq data=dtemp.sim_student_flag;
tables flag_resid;
run;

*15% is high...;
ods pdf close ; 

/* save these frquencies in documentation but
only pvalue is needed for next step */

ods exclude all;

proc ttest data= dTemp.sim ; 
by featureid;
class sex;
var log_apn_uq_ff;
ods output  equality=dtemp.sim_equality ttests=dtemp.sim_ttests ;
run;

ods exclude none;
ods results;

