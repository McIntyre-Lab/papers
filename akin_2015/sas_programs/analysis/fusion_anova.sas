/* import libraries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* subset for ANOVAs */

data counts_all_on;
  set sugrue.counts_w_means;
  if flag_control_on=1 and flag_treat_on=1 then output;
run;

/* Run ANOVA */

* sort data;
proc sort data= counts_all_on  nodup;
  by fusion_id subject;
  run;

/* freq for order */
proc freq data= counts_all_on;
   tables treatment;
run;

/* ANOVA - took a while to run (~17hr) not sure why */

ods listing close;
	ods pdf close;
	ods html close;
	ods graphics off;
proc glm data= counts_all_on;
by fusion_id;
class treatment;
model log_apn=treatment;  

lsmeans treatment;			/*0	1*/
contrast '0-1' treatment		-1	1;
estimate '0-1' treatment		-1	1;

ods output modelanova = anova
	lsmeans = lsmean
	contrasts = contrasts
	estimates = estimates;
	output out=resid predicted=pred residual=resid student=student;
run;
quit;
ods listing;


/* Flag residuals that fail normality */
ods graphics off;
proc univariate data = resid normal noprint;
  by fusion_id;
  var Resid;
  output out = normtest probn=pnorm;
  run;


data sugrue.flag_resids_logapn;
  set normtest;
  if pnorm ge 0 and pnorm le 0.05 then flag_fail_norm = 1;
        else if pnorm='.' then flag_fail_norm=.;
	else flag_fail_norm = 0;
  run;

ods graphics on;
ods listing;
	ods pdf;
	ods html;

proc freq data = sugrue.flag_resids_logapn noprint;
  tables flag_fail_norm / out=resid_fail_norm ;
  run;


/* flag fail norm summary 

0 = 182572 (96.08%)
1 = 7449 (3.92%)

*/ 


/* ANOVA */


ods listing close;
	ods pdf close;
	ods html close;
	ods graphics off;
proc glm data=counts_all_on;
by fusion_id;
class treatment;
model log_apn= treatment;  
lsmeans treatment;	
ods output modelanova = anova
	lsmeans = lsmean;
	output out=resid predicted=pred residual=resid student=student;
run;
quit;
ods listing;

data sugrue.anova_logapn ;
set anova ;
run ;
data sugrue.lsmean_logapn ;
set lsmean ;
run ;



/*grab only on fusions */

data fusions_all_on;
 set counts_all_on;
 keep fusion_id;
run;

proc sort data=fusions_all_on nodup;
by fusion_id;
run;

proc sort data=sugrue.flag_resids_logapn;
   by fusion_id;
run;

proc sort data=sugrue.anova_logapn;
   by fusion_id;
run;


data sugrue.flag_resids_2;
   merge sugrue.flag_resids_logapn (in=in1) fusions_all_on (in=in2);
   by fusion_id;
   if in1 and in2 then output;
run;

data sugrue.anova_2;
   merge sugrue.anova_logapn (in=in1) fusions_all_on (in=in2);
   by fusion_id;
   if in1 and in2 then output;
run;
