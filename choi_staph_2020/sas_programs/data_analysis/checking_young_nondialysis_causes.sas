
proc means data=choi.choi_756 min p25 p50 p75 max;
where sab_status="R" and dialdep=1;
class race_dialdep;
var age;
run;



proc means data=model min p25 p50 p75 max;
class  sab_status dialdep;
var age;
run;

proc freq data=model;
where sab_status="R" and dialysis_race="zilch_dial";
tables race_2level;
run;


proc glm data=model;
class dialysis_race sab_status;
model age=dialysis_race|sab_status;
run;


proc freq data=model;
where sab_status="R" and dialdep=0;
tables foreignb surgery30d;
run;


proc freq data=model;
where sab_status="R" and dialdep=0;
tables flag_age_above_med*(
neoplasm                
transpat                
steroid                 
hiv                     
foreignb   
route 
surg30d            
Persistent   
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept)/exact;  
run;


data peeking;
set model;
if mrsa=1 or metendo=1 or foreignb=1 then common=1;
else common=0;

if common=1 then metveggie="mvfb";
else metveggie="nope";
run;

proc freq data=peeking;
where sab_status="R" ;
tables flag_age_above_med*common*dialysis_race;
run;


proc freq data=peeking;
where sab_status="R" and flag_age_above_med=0;
tables common*race_2level;
run;

data subset;
set model;
where sab_status="R" and dialdep=0
and race=3;
run;

proc freq data=subset;

tables 
flag_age_above_med
neoplasm                
transpat                
steroid                 
hiv                     
foreignb   
route 
surg30d            
Persistent   
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept;  
run;

proc freq data=peeking;
where flag_age_above_med=1 and common=1;
tables race_2level*outcome;
run;


proc freq data=peeking;
tables flag_age_above_med*common;
run;


proc freq data=peeking;
tables common*outcome;
run;


proc logistic data=peeking;
class dialysis_race flag_apache_above_med common metveggie;
model flag_rsab=dialysis_race 
		 flag_apache_above_med metveggie /clodds=pl;
 ods output fitstatistics=fit_red modelanova=model_red cloddspl=clodds_red;

		run;
