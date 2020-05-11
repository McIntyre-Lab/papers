libname relapse "C:\a1stuff\staph_relapse";

options ls=110 ps=54 errors=2 cleanup nofmterr pageno=1 mergenoby=warn noovp;

proc format library=relapse;
   value ynuf    0 = 'No'
                 1 = 'Yes';
   ;
   value genderf 0 = 'Male'
                 1 = 'Female'
   ;
   value routef  1 = 'Hospital'
                 3 = 'Healthcare'
                 4 = 'Community'
   ;
   value racef   1 = 'White'
                 2 = 'Black'
                 3 = 'Other'
   ;
   value outcof  1 = 'Alive'
                 2 = 'Recur Infection'
                 3 = 'Died - Gram neg'
                 4 = 'Died - Other'
   ;
run;
 
data relapse;
set relapse.choi_19;
run;

proc contents data=relapse;
run;

proc freq data=relapse;
tables choi_era;
run;

proc freq data=relapse;
*where choi_era=1;
*tables patient_number/out=count_patient;
tables choi_unique choi_era choi;
run;


proc freq data=relapse;
where choi=1;

tables patient_number/out=count_patient;
run;


proc univariate data=relapse;
 var age;
 run;

 data check_3838;
 set relapse; where id=3838;
 run;

 proc freq data=relapse;
 tables pos_cult_dt;
 run;


 proc freq data=relapse ;
 *tables patientnumber/out=count_patients;
 tables site;
 run;

 proc freq data=relapse;
 where site ne "";
 tables choi_era;
 run;

 data find_ssab;
 set relapse;
 where choi_era=1;
 run;


proc freq data=find_ssab;
*tables outcome;
*tables route;
tables recur_sa;
run;


proc freq data=relapse;
*tables outcome;
*tables route;
tables recur_sa*choi_era;
run;

data s_sab_clin;
set relapse;
where choi_era=1 and recur_sa=0;
Sab_status="S";
if outcome=. then delete;
run;

proc sort data=relapse;
by id;

proc sort data=pfge;
by id;
run;

data r_sab_clin not_r oops;
merge relapse (in=in1) pfge(in=in2);
by id;
if in1 and in2 then output r_sab_clin;
else if in1 then output not_r;
else if in2 then output oops;
run;

data choi;
set s_sab_clin r_sab_clin;
if sab_status="" then sab_status="R";
run;

proc freq data=choi;
tables sab_status;
run;

data check_5588;
set relapse;
where id=5588;
run;

proc freq data= s_sab_clin;
tables outcome;
run;


proc sort data=relapse;
by id;
run;

data check;
set relapse;
if id >4600 and id <4800;
run;

