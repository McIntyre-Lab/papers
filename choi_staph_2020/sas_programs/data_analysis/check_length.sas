

proc sort data=relapse;
by id;

proc sort data=pfge_included;
by id;

data clin_pfge;
merge relapse (in=in1) pfge_included (in=in2);
by id;
if in2;
run;

proc freq data=clin_pfge;
tables patientnumber/out=count_rsab;
run;


/*69 patients!*/


proc univariate data=clin_pfge normal plot;
var days;
run;

proc sort data=clin_pfge;
by pfge;

proc univariate data=clin_pfge normal plot;
where pfge ne "";
by pfge;
var days;
run;

proc ttest data=clin_pfge;
class pfge;
var days;
run;

proc npar1way data=clin_pfge;
class pfge;
var days;
run;

/*by patient*/

proc sort data=clin_pfge;
by patientnumber;

proc sort data=patient_pfge;
by patientnumber;

data patient_classify;
merge clin_pfge patient_pfge;
by patientnumber;
if days=. then delete;
run;

/* 85 pairs*/

proc univariate data=patient_classify;
where num_pairs=1;
var days;
run;


proc freq data=patient_classify;
where num_pairs=1 and days le 150;
tables pfge;
run;

proc freq data=patient_classify;
where num_pairs=1 and days le 60;
tables pfge;
run;




