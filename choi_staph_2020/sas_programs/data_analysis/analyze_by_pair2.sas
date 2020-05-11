
/*ctokine paired analysis*/
/* */


data choi.diff_pair;
set choi.set_pair;
diff_rantes= recur_rantes-control_rantes;
reinfect_early=reinfect;
if days le 60 then reinfect_early=2;
run;


proc freq data=diff_pair;
tables pfge num_pairs;
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
tables recur_outcome*recur_rantes_class;
run;


proc export data=diff_pair
            OUTFILE="Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis\diff_pair.csv"

            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
