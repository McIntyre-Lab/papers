libname relapse "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\relapse_data_from_duke";
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";
libname wgs "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata";

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
 /*formats error so clear formats*/

proc datasets lib=relapse memtype=data;
   modify choi_19;
     attrib _all_ label=' ';
     attrib _all_ format=;
contents data=relapse.choi_19;
run;
quit;

proc contents data=relapse.choi_19;
run;

/*note that the spatype flags in the original data are inconsistent*/
data relapse;
set relapse.choi_19;
if ridom_cc="2" then flag_ridom_cc2=1;
	else if ridom_cc="" then flag_ridom_cc2=.;
	else flag_ridom_cc2=0;
if ridom_cc="8" then flag_ridom_cc8=1;
	else if ridom_cc="" then flag_ridom_cc8=.;
	else flag_ridom_cc8=0;
if ridom_cc="2" then ridom_cc2_vs_8=2;
	else if ridom_cc="8" then ridom_cc2_vs_8=8;
	else if ridom_cc="" then ridom_cc2_vs_8=.;
	else ridom_cc2_vs_8=0;

	if apache_total ge 16 then flag_apache_above_med=1;
	else if apache_total=. then flag_apache_above_med=.;
	else flag_apache_above_med=0;
if age ge 60  then flag_age_above_med=1;
	else if age=. then flag_age_above_med=.;
	else flag_age_above_med=0;
if outcome=1 then flag_alive=1;
else if outcome=. then flag_alive=.;
else if outcome=2 then flag_alive=1;
else flag_alive=0;

if race=3 and dialdep=1 then dialysis_race="AA_dialdep";
	else if race=2 and dialdep=1 then dialysis_race="C_dialdep";
	else if race=. or dialdep=. then dialysis_race="";
	else if dialdep=1 then dialysis_race="O_dialdep";
	else dialysis_race="zilch_dialysis";

/*if surg30d=1 and foreignb=1 then surg30d_fb="surg_and_fb";
	else this was not helpful to separate*/ 
	if surg30d=1 or foreignb=1 then surg30d_fb="surg_OR_fb";
	else if surg30d=. or foreignb=. then surg30d_fb="";
	else surg30d_fb="surg_zilch_fb";

	if  flag_age_above_med=1 and flag_apache_above_med=1 then apache_age="sick_and_old";
		else if  flag_age_above_med=1 or flag_apache_above_med=1 then apache_age="sick_or_old";
		else if flag_age_above_med=. or flag_apache_above_med=. then apache_age="";
		else apache_age="zilch";

keep id age apache_total 
female race 
dm  dialdep idu neoplasm                
transpat steroid hiv foreignb                
surg30d Persistent mrsa metendo metabsc metarth metsept
outcome choi_era recur_sa
route
RIDOM_CC flag_ridom_cc2 flag_ridom_cc8 ridom_cc2_vs_8
LOS flag_alive 
dialysis_race surg30d_fb flag_age_above_med flag_apache_above_med
;
run;

data s_sab_clin_686;
set relapse;
where choi_era=1 and recur_sa=0; /*see the look programthis is best guess*/
if outcome=. then delete;
run;

data isolate_5588;
set relapse;
where id=5588;
outcome=1; /*updated because of note from felicia ruffin*/
run;

/* 3/29/2020 Hey Lauren,

Correct, 5588 was a single episode. He was coded as recurrence due to a single blood culture within the 90-day window. It was not positive Staph aureus (SA), but Staph coag neg (SCN).  The positive culture was likely a contaminant because the SCN was in one culture only. The outcome should be changed to "Cure".  If it had been in both culture, it would be an outcome of "recurrence" and the next question asks if it was Staph aureus which would have been answered as "no".  It is kind of confusing because that outcome variable encompasses a little more than the way "recurrence" is define in the study overall.  
*/

data s_sab_clin_687;
set s_sab_clin_686 isolate_5588;
Sab_status="S";
run;

proc sort data=relapse;
by id;
/*from program analyze_pfge.sas*/
proc sort data=choi.id_link_2_patient_pfge;
by id;

data r_sab_clin not_r oops;
merge relapse (in=in1) choi.id_link_2_patient_pfge(in=in2);
by id;
if in1 and in2 then output r_sab_clin;
else if in1 then output not_r;
else if in2 then output oops;
run;

data choi;
set s_sab_clin_687 r_sab_clin;
if sab_status="" then sab_status="R";
run;

data wgs;
set wgs.sampleid_st_cc_ref_kitchen_sink;
rename sampleid=id;
drop pairnum pfge days episode_number episodes count;
run;

proc sort data=wgs;
by id;
proc sort data=choi;
by id;
proc sort data=cytokines;
by id;

data choi_wgs;
merge choi(in=in1) wgs (in=in2) cytokines(in=in3);
by id;
run;


data choi.choi_841;
set choi_wgs;
if patientnumber=0 and id=4132 then episode_number=1;
if patientnumber=0 and id=4130 then episode_number=2;

run;

data rsab_nopfge;
set r_sab_clin;
drop pfge days episodes num_pairs;
if patientnumber=0 and id=4132 then episode_number=1;
if patientnumber=0 and id=4130 then episode_number=2;
run;

proc sort data=rsab_nopfge out=r_sab_clin_sort ;
by patientnumber episode_number;
run;

data r_sab_69_unique;
set r_sab_clin_sort ;
by patientnumber;
if first.patientnumber;
run;


/*sa type seems randomly assigned to isolates */

data choi_756_clin;
set s_sab_clin_687 r_sab_69_unique;
if sab_status="" then sab_status="R";

if race=2 then race_3level="Caucasian";
else if race=3 then race_3level="African_American";
else if race=. then race_3level="";
else race_3level="Other";

if race=2 then race_2level="Caucasian";
else if race=3 then race_2level="African_American";
else race_2level="";
run;

proc sort data=choi_756_clin;
by patientnumber;
proc sort data=choi.patient_pfge_all_episodes;
by patientnumber;

data choi_756;
merge  choi_756_clin(in=in1) choi.patient_pfge_all_episodes (in=in2);
by patientnumber;
run;

data choi.choi_756;
set choi_756;

run;


proc export data=choi.choi_841
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis\choi_841.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


proc export data=choi.choi_756
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis\choi_756.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

