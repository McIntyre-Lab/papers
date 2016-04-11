/* Proc glm and proc mixed to look at residuals..*/

libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata' ;
libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data' ;

proc sort data=ribo.all_coverage_counts_with_key ;
by fusion_id;
run;
proc sort data=ribo.on_calls_gt_apn0;
by fusion_id;
run;



data ribo_all ;
merge ribo.all_coverage_counts_with_key (in=in1) ribo.on_calls_gt_apn0;
by fusion_id ;
run;

/*keep fusions that are on for every treatment*/
data ribo_allfusions_alltrt_on;
set ribo_all ;
if flag_fusion_all_on0 = 1 ;
run;

proc freq data = ribo_allfusions_alltrt_on;
tables flag_fusion_all_on0;
run;

/*proc glm on treatment, output residuals*/

proc sort data = ribo_allfusions_alltrt_on;
by fusion_id;
run;

ods listing close ;

proc glm data = ribo_allfusions_alltrt_on ;
by fusion_id ;
class trt ;
model logRPKM = trt ;
output out= resid r = resid p=pred ;
ods output modelanova = model ;
run ;
quit ;

ods listing;

/*univariate test to look at residuals. Residuals should be less than 10% of total that fail normality*/
proc univariate data = resid normal noprint ;
by fusion_id ;
var resid ;
output out = normtest probn=pnorm ;
run;


data normflags ;
set normtest ;
if pnorm = . then flag_fail_norm = .;
	else if pnorm le 0.05 then flag_fail_norm = 1 ;
	else flag_fail_norm = 0;
run ;


ods listing close;
/*keep in mind that this is different from glm this is the mixed model*/
proc mixed data=ribo_allfusions_alltrt_on;
by fusion_id ;
class trt;
model logRPKM= trt / htype=1 outp=mixresid;
lsmeans trt/diff ;

ods output tests1 = mixedon
           lsmeans = mixedlsmean
           diffs = mixeddiffs
;
run; quit;
/*make sure you open again*/
ods listing;


/*look at residuals from the mixed model */
proc univariate data = mixresid normal noprint;
by fusion_id;
var Resid;
output out = normtest probn = pnorm;
run;

data mix_normflags ;
set normtest;
