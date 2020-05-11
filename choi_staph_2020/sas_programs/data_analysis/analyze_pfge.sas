
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";

data pfge_included;
set choi.pfge_included;

proc freq data=pfge_included;
tables patientnumber*PFGE/out=count_pfge;
run;

/*54 pairs identical, 31 pairs differnt*/

proc freq data=count_pfge;
tables  patientnumber/out=count_pfge_results;
run;

data list_multiple_pfge_results; /*note these are all multiple episodes, this is just where the PFGE results change*/ 
set count_pfge_results;
where COUNT >2;
run;

/*6 multiple this means that there are differnces among the isolate types for multiple episodea*/
proc sort data=list_multiple_pfge_results;
by patientnumber;

proc sort data=pfge_included;
by patientnumber;

data multiple single;
merge pfge_included (in=in1) list_multiple_pfge_results(in=in2);
by patientnumber;
if in2 then output multiple;
else output single;
run;

/*135 episodes in single; 19 in multiple*/

proc freq data=single;
where pfge ne "";
tables patientnumber*PFGE/out=count_single_pfge;
run;

/*47 identical pairs and 25 differnt pairs)*/

proc freq data=count_single_pfge;
tables pfge*count;
run;

/* 63- patients : 55 have 2 episodes (1 pair), 7 have 3 episodes (2 pairs), 1 have 4 episodes (3 pairs)
totzl is 72 pairs*/


proc freq data=multiple;
where pfge ne "";
tables patientnumber*PFGE/out=count_multiple_pfge;
run;

/*13 pairs 6 differnt 7 identical*/
/*6 patients: 5 with 2 pairs one identical and 1 differnt;
1 with 3 pairs 2 identical 1 different*/


/*  47+7 identical pairs , 25+6 differnt pairs total is  85 pairs*/



/*classify by patient*/

data single_describe;
set count_single_pfge;
rename count=num_pairs pfge=episodes;
drop percent;
run;

proc sort data=count_multiple_pfge;
by patientnumber;

/* for the mutiple types make sure pair count is accurate by setting manually*/
data multiple_describe;
set count_multiple_pfge;
by patientnumber;
if first.patientnumber;
episodes="Mixed";
if patientnumber=65 then num_pairs=3; else num_pairs=2;
drop pfge percent count;

run;

data patient_pfge;
set single_describe multiple_describe;
run;

/*check 85 pairs- yes!*/
proc means data=patient_pfge sum;
var num_pairs;
run;

proc sort data=pfge_included;
by patientnumber id;

data pfge_include_epnum;
set pfge_included;
episode_number + 1;
by patientnumber;
if first.patientnumber then episode_number = 1;
run ;


proc sort data=patient_pfge;
by patientnumber ;

proc freq data=patient_pfge;
tables num_pairs;
run;

/* this is the stack file with the isolate id linked to a patient number, 
number of pairs and the classification of the patient based on the set of pairs*/

data id_link_2_patient_pfge oops1 oops2;
merge pfge_include_epnum (in=in1) patient_pfge (in=in2);
by patientnumber ;
if in1 and in2 then output id_link_2_patient_pfge;
else if in1 then output oops1;
else if in2 then output oops2;
run;

data choi.id_link_2_patient_pfge;
set id_link_2_patient_pfge;
run;

proc export data=id_link_2_patient_pfge
            OUTFILE= "C:\a1stuff\staph_relapse\pfge_results.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


proc sort data=id_link_2_patient_pfge;
by patientnumber id;
run;


data one_pair two_or_more_pairs;
set id_link_2_patient_pfge;
if days=. then delete;
if num_pairs=1 then output one_pair;
else output two_or_more_pairs;
run;

data two_or_more_pairs1;
set  two_or_more_pairs;
keep patientnumber episodes num_pairs;

run;

/*14 patients with 2 or more pairs*/
proc sort data=two_or_more_pairs1 out=two_or_more_pairs_classification nodupkey;
by patientnumber episodes num_pairs;
run;


proc transpose data=two_or_more_pairs out=days_flip;
by patientnumber;
var days;
run;

data days_by_episode;
set days_flip;
rename col1=days 
col2=days2
col3=days3;
drop _name_;
run;

proc transpose data=two_or_more_pairs out=pfge_flip;
by patientnumber;
var pfge;
run;

data pfge_by_episode;
set pfge_flip;
rename col1=pfge 
col2=pfge2
col3=pfge3;
drop _name_;
run;

data patient_multi_episode;
merge two_or_more_pairs_classification days_by_episode pfge_by_episode;
by patientnumber;
run;

data choi.patient_pfge_all_episodes;
set patient_multi_episode one_pair;
drop id episode_number;

if pfge="Different" then reinfect=1;
else if days>150 then reinfect=1;
	else reinfect=0;


if pfge="Different" then reinfect_relapse=2;
else if days>150 then reinfect_relapse=2;
	else if days le 60 then reinfect_relapse=1;
	else reinfect_relapse=0;

if pfge="Different" then early_relapse=0;
else if days>150 then early_relapse=0;
	else if days le 60 then early_relapse=1;
	else early_relapse=0;
run;


proc export data=choi.patient_pfge_all_episodes
            OUTFILE= "C:\a1stuff\staph_relapse\patient_pfge_all_episodes" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc freq data=patient_pfge_all_episodes;
tables pfge*reinfect;
run;







