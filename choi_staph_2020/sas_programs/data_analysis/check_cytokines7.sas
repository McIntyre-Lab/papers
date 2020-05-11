
proc sort data=cytokines;
by id;
proc sort data=choi.choi_756;
by id;


data cyto_clin no_cyto oops;
merge choi.choi_756 (in=in1) cytokines(in=in2);
by id;
if in1 and in2 then output cyto_clin; /*should be and is 42*/
else if in1 then output no_cyto;
else  output oops;
run;

proc freq data=cyto_clin;
tables gender female recur_sa recurren recurrent;
run;

/* very strange that recur_sa recurren recurrent do not gives the right (21) count*/
/*gender and female both give 19 females and 23 males meaning these are not matched completely on sex*/



/*
proc sort data=pfge_included;
by id;
*/
proc sort data=choi.id_link_2_patient_pfge;
by id;

data id;
set choi.id_link_2_patient_pfge;
keep id;
run;
/*note pfge info now is clin data!*/

data cyto_clin_pfge cyto_control pfge_only;
merge cyto_clin (in=in1) id (in=in2);
by id;
if in1 and in2 then output  cyto_clin_pfge; /*21*/ 
else if in1 then output cyto_control; /*21*/ 
else if in2 then output pfge_only; /*133*/
run;

/*figure out the match between the recur and recurrent*/

data cyto_control1;
set cyto_control;
control_rantes=rantes;
control_id=id;
control_outcome=outcome;
control_dialdep=dialdep;
control_race=race;
control_age=age;
control_female=female;
control_apache_total=apache_total;
control_ridom_cc=ridom_cc;
control_mrsa=mrsa;
if RIDOM_CC="" then RIDOM_CC="NT";
controlnum=_n_;
keep age female race control_rantes RIDOM_CC controlnum 
control_id control_outcome control_race control_dialdep control_apache_total
control_age control_female control_ridom_CC control_age control_female control_mrsa;
run;

data cyto_clin_pfge1;
set cyto_clin_pfge;
recur_rantes=rantes;
rantes_id=id;
recur_outcome=outcome;
recur_dialdep=dialdep;
recur_race=race;
recur_age=age;
recur_female=female;
recur_apache_total=apache_total;
recur_ridom_cc=ridom_cc;
recur_mrsa=mrsa;
if RIDOM_CC="" then RIDOM_CC="NT";
keep age female race recur_rantes rantes_id 
patientnumber days days2 days3 pfge pfge2 pfge3 
episodes num_pairs reinfect reinfect_early
RIDOM_CC recur_outcome recur_race 
recur_dialdep recur_apache_total
recur_age recur_female recur_ridom_CC recur_mrsa;
run;
/*find reinfect_early*/

proc freq data=cyto_clin_pfge1;
tables num_pairs;
run;

proc sort data=cyto_control1;
by  RIDOM_CC race  female age;
run;

proc sort data=cyto_clin_pfge1;
by   RIDOM_CC race  female age;
run;

data pair_rantes nomatch_control nomatch_case;
merge cyto_control1 (in=in1) cyto_clin_pfge1 (in=in2);
by RIDOM_CC race  female;
if in1 and in2 then output pair_rantes;
else if in1 then output nomatch_control;
else if in2 then output nomatch_case;

run;

/* patientnum 82 matchs both 5 and 13; 5 is the same age;*/

data pair_rantes2;
set pair_rantes;
if controlnum=13 then delete;
drop RIDOM_CC race  female age;
run;

data mismatch;
set  cyto_control1 ;
if controlnum=13;

/*look by eye at the 3 mismatches*/
/*match on race and age*/
/*in order by coincidence*/
data controls_left;
set mismatch nomatch_control;
match=_n_;

keep match control_rantes controlnum control_id control_outcome;
run;

data case_left;
set nomatch_case;
match=_n_;

drop RIDOM_CC race  female age;
run;
proc sort data=controls_left ;
by match;
proc sort data=case_left;
by match;

data mismatch3;
merge   case_left (in=in2) controls_left (in=in1);
by match;
if in1 and in2;
run;

/*put pairs back together*/

data choi.set_pair;
set pair_rantes2 mismatch3;
run;


proc means data=set_pair median;
var control_rantes recur_rantes;
run;


/*check the median before the split this should not change*/
data cyto_control1;
set cyto_control;
sab_status="S";
keep  id sab_status rantes;
run;

data cyto_clin_pfge1;
set cyto_clin_pfge;
sab_status="R";
keep  id sab_status rantes;
run;

/* need to fix

data cytokines;
set cyto_control1 cyto_clin_pfge1;
drop control_: recur_: ;
run;

proc sort data=choi.choi_756;
by id;

proc sort data=cytokines;
by id;

data clin_cyto_check_outcome;
merge choi.choi_756 (in=in1) cytokines(in=in2);
by id;
if in1 and in2;
run;

proc sort data=clin_cyto_check_outcome;
by outcome;
proc univariate data=clin_cyto_check_outcome normal plot;
by outcome;
var rantes;
run;

proc freq data=clin_cyto_check_outcome;
tables sab_status*outcome;
run;

proc glm data=clin_cyto_check_outcome;
class sab_status;
model rantes=sab_status;
run;

proc glm data=clin_cyto_check_outcome;
where outcome=1;
class sab_status;
model rantes=sab_status;
run;

proc npar1way data=clin_cyto_check_outcome;
where outcome=1;
class sab_status;
var rantes;
run;


proc print data=clin_cyto_check_outcome;
where sab_status="S";
var id outcome;
run;

proc  npar1way data=cytokines;
class sab_status;
var rantes;
run;
