
/*ctokine paired analysis*/
/* */


data diff_pair;
set set_pair;
diff_rantes= recur_rantes-control_rantes;
diff_EOTAXIN	=	recur_EOTAXIN	-	control_EOTAXIN	;
diff_IFNa	=	recur_IFNa	-	control_IFNa	;
diff_IFNy	=	recur_IFNy	-	control_IFNy	;
diff_IL6	=	recur_IL6	-	control_IL6	;
diff_IL8	=	recur_IL8	-	control_IL8	;
diff_IL12	=	recur_IL12	-	control_IL12	;
diff_IL13	=	recur_IL13	-	control_IL13	;
diff_IL15	=	recur_IL15	-	control_IL15	;
diff_IL1B	=	recur_IL1B	-	control_IL1B	;
diff_IL1RA	=	recur_IL1RA	-	control_IL1RA	;
diff_IL2R	=	recur_IL2R	-	control_IL2R	;
diff_IP10	=	recur_IP10	-	control_IP10	;
diff_MCP1	=	recur_MCP1	-	control_MCP1	;
diff_MIG	=	recur_MIG	-	control_MIG	;
diff_MIP1B	=	recur_MIP1B	-	control_MIP1B	;

if	diff_EOTAXIN	>0	then 	posdiff_EOTAXIN	=	1;	else	posdiff_EOTAXIN	=	0;
if	diff_IFNa	>0	then 	posdiff_IFNa	=	1;	else	posdiff_IFNa	=	0;
if	diff_IFNy	>0	then 	posdiff_IFNy	=	1;	else	posdiff_IFNy	=	0;
if	diff_IL6	>0	then 	posdiff_IL6	=	1;	else	posdiff_IL6	=	0;
if	diff_IL8	>0	then 	posdiff_IL8	=	1;	else	posdiff_IL8	=	0;
if	diff_IL12	>0	then 	posdiff_IL12	=	1;	else	posdiff_IL12	=	0;
if	diff_IL13	>0	then 	posdiff_IL13	=	1;	else	posdiff_IL13	=	0;
if	diff_IL15	>0	then 	posdiff_IL15	=	1;	else	posdiff_IL15	=	0;
if	diff_IL1B	>0	then 	posdiff_IL1B	=	1;	else	posdiff_IL1B	=	0;
if	diff_IL1RA	>0	then 	posdiff_IL1RA	=	1;	else	posdiff_IL1RA	=	0;
if	diff_IL2R	>0	then 	posdiff_IL2R	=	1;	else	posdiff_IL2R	=	0;
if	diff_IP10	>0	then 	posdiff_IP10	=	1;	else	posdiff_IP10	=	0;
if	diff_MCP1	>0	then 	posdiff_MCP1	=	1;	else	posdiff_MCP1	=	0;
if	diff_MIG	>0	then 	posdiff_MIG	=	1;	else	posdiff_MIG	=	0;
if	diff_MIP1B	>0	then 	posdiff_MIP1B	=	1;	else	posdiff_MIP1B	=	0;
if	diff_rantes	>0	then 	posdiff_rantes	=	1;	else	posdiff_rantes	=	0;

reinfect_early=reinfect;
if days le 60 then reinfect_early=2;
run;


proc freq data=diff_pair;

tables pfge num_pairs;
tables posdiff_EOTAXIN
posdiff_IFNa
posdiff_IFNy
posdiff_IL6
posdiff_IL8
posdiff_IL12
posdiff_IL13
posdiff_IL15
posdiff_IL1B
posdiff_IL1RA
posdiff_IL2R
posdiff_IP10
posdiff_MCP1
posdiff_MIG
posdiff_MIP1B
posdiff_rantes;

run;

/*overall test from the paper*/


title1 "control vs Rsab";
proc univariate data=diff_pair;
var diff_rantes;
run;
/*signed rank 0.0053*/
*checking influence of the 3 going negative;


title1 "control vs Rsab no neg diff";
proc univariate data=diff_pair;
where diff_rantes>0;
var diff_rantes;
run;
/*signed rank <0.0001*/

title1 "full set n=21";
proc means data=diff_pair p25 p50 p95;
var control_rantes recur_rantes diff_rantes;
run;

proc means data=diff_pair p25 p50 p95;
class reinfect;
var control_rantes recur_rantes diff_rantes;
run;

title1 "reinfect_early =2 is recurr before 60 days";
proc means data=diff_pair p25 p50 p95;
class reinfect_early;
var control_rantes recur_rantes diff_rantes;
run;

title "all pfge=different are more than 60 days";
proc means data=diff_pair p25 p50 p95;
class pfge;
var control_rantes recur_rantes diff_rantes;
run;



/* checking for sample not in original 686 is not influential all looks 0k*/
/*proc univariate data=diff_pair;
where control_id ne 5588;
var diff_rantes;
run;

proc ttest data=diff_pair;
paired control_rantes*recur_rantes;
run;

proc ttest data=diff_pair;
where control_id ne 5588;
paired control_rantes*recur_rantes;
run;*/


/*split by recur, reinfect early*/

title1 "recurr vs reinfect";

proc sort data=diff_pair;
by reinfect;

proc npar1way data=diff_pair;
class reinfect;
var diff_rantes;
run;

proc means data=diff_pair p25 p50 p75;
class reinfect;
var control_rantes recur_rantes diff_rantes;
run;

title1 "recurr vs reinfect vs recurr_early";
proc sort data=diff_pair;
by reinfect_early;


proc means data=diff_pair p25 p50 p75;
class reinfect_early;
var control_rantes recur_rantes diff_rantes;
run;



proc means data=diff_pair p25 p50 p75;
class control_outcome;
var control_rantes recur_rantes diff_rantes;
run;

title1 "test pfge status";
proc sort data=diff_pair;
by pfge;

proc univariate data=diff_pair normal plot;
by pfge;
var diff_rantes;
run;

proc freq data=diff_pair;
tables pfge*class_recur;
run;

proc sort data=diff_pair;
by class_recur;

proc ttest data=diff_pair;
by pfge;
paired control_rantes*recur_rantes;
run;

proc univariate data=diff_pair normal plot;
by class_recur;
var diff_rantes;
run;

proc gplot data=diff_pair;
plot days*diff_rantes=control_outcome;
plot days*diff_rantes=class_recur;
plot days*recur_rantes;
run;


proc sort data=diff_pair;
by recur_rantes_high;
proc univariate data=diff_pair normal plot;
by recur_rantes_high;
var days;
run;


proc sort data=diff_pair;
by recur_outcome;
proc univariate data=diff_pair normal plot;
by recur_outcome;
var days;
run;



proc sort data=diff_pair;
by pfge;
proc univariate data=diff_pair normal plot;
where control_id ne 5588;
by pfge;
var days;
run;

proc glm data=diff_pair;
where pfge="Identical";
model days=recur_rantes;
run;

proc glm data=diff_pair;
where pfge="Different";
model days=recur_rantes;
run;

proc glm data=diff_pair;

model days=recur_rantes;
run;


proc npar1way data=diff_pair;
where control_id ne 5588;
class pfge;
var recur_rantes;
run;

proc glm data=diff_pair;

class recur_mult;
model recur_rantes=recur_mult;
run;


proc glm data=diff_pair;
class recur_rantes_class;
model days=recur_rantes_class;
run;

proc means data=diff_pair median;
class recur_mult;
var recur_rantes;
run;

proc freq data=diff_pair;
tables recur_outcome*recurent;
run;


proc export data=diff_pair

            OUTFILE="Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis\diff_pair.csv"

            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
