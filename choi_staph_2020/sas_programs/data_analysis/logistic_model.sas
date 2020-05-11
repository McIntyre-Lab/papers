
/*
apache_total 
age*/
proc means data=choi_756_clin med;
var apache_total age;
run;

data model;
set choi_756_clin;

where race=2 or race=3;
if apache_total ge 16 then flag_apache_above_med=1;
	else if apache_total=. then flag_apache_above_med=.;
	else flag_apache_above_med=0;
if age ge 60  then flag_age_above_med=1;
	else if age=. then flag_age_above_med=.;
	else flag_age_above_med=0;
if sab_status="R" the flag_rsab=0;
	else flag_rsab=1;
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
run;

proc freq data=model;
tables sab_status*(flag_apache_above_med flag_age_above_med)/exact;
run;

/*check for confounding*/
proc freq data=model;
tables dialdep*(race_2level foreignb surg30d persistent neoplasm steroid flag_apache_above_med flag_age_above_med)/exact;
tables race_2level *(foreignb surg30d persistent neoplasm steroid flag_apache_above_med flag_age_above_med)/exact;
tables foreignb*(surg30d persistent neoplasm steroid flag_apache_above_med flag_age_above_med)/exact;
tables surg30d *(persistent neoplasm steroid flag_apache_above_med flag_age_above_med)/exact;
tables persistent*( neoplasm steroid flag_apache_above_med flag_age_above_med)/exact;
tables neoplasm *(steroid flag_apache_above_med flag_age_above_med)/exact;
tables steroid* (flag_apache_above_med flag_age_above_med)/exact;
tables flag_apache_above_med*flag_age_above_med/exact;
run;

/*unadjusted effects*/


%unadjusted (flag_apache_above_med);

data fit_unadjusted;
set fit_flag_apache_above_med;


data model_unadjusted;
set model_flag_apache_above_med;


data clodds_unadjusted;
set clodds_flag_apache_above_med;
run;


%macro unadjusted (effect);
proc logistic data=model;
class dialdep  race_2level foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med mrsa 
	dialysis_race  surg30d_fb apache_age;
	model flag_rsab=&effect/clodds=PL;
	ods output fitstatistics=fit_&effect  modelanova=model_&effect cloddspl=clodds_&effect;
		run;
	
data fit_unadjusted;
	set fit_unadjusted fit_&effect;
	name="&effect";
	run;

data model_unadjusted;
	set model_unadjusted model_&effect;
	run;

data clodds_unadjusted;
	set clodds_unadjusted clodds_&effect;
	run; 
%mend;

%unadjusted (dialdep);
%unadjusted (race_2level);
%unadjusted (foreignb);
%unadjusted (surg30d); 
%unadjusted (persistent);
%unadjusted (neoplasm);
%unadjusted (steroid);
%unadjusted (mrsa);
%unadjusted (flag_age_above_med);
%unadjusted (dialysis_race);

/*%unadjusted (surg30d_fb); diveds the data into old and sick vs surg/fb*/
/*%unadjusted(apache_age); strange combo*/

   
proc freq data=model;
tables dialdep*race  flag_age_above_med*flag_apache_above_med ;
run;

proc export data=model_unadjusted
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\logistic\model_unadjusted.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc export data=clodds_unadjusted
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\logistic\clodds_unadjusted.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


		/*full model*/


proc logistic data=model;
class   dialdep  race_2level foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med mrsa 
	dialysis_race  surg30d_fb apache_age;
model flag_rsab=persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med mrsa 
	dialysis_race   foreignb surg30d /clodds=PL;
ods output fitstatistics=fit_full modelanova=model_full cloddspl=clodds_full;
run;



proc export data=model_full
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\logistic\model_full.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc export data=clodds_full
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\logistic\clodds_full.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


/*drop steroid*/
proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race foreignb surg30d persistent neoplasm 
		 flag_apache_above_med flag_age_above_med mrsa;
		run;
/*drop neoplasm*/
		
proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race foreignb surg30d persistent 
		 flag_apache_above_med flag_age_above_med mrsa;
		run;

/*drop persistent*/

		
proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race foreignb surg30d 
		 flag_apache_above_med flag_age_above_med mrsa;
		run;

/*drop age*/

proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race foreignb surg30d 
		 flag_apache_above_med mrsa /clodds=pl;
		run;

		/*drop mrsa*/
proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race foreignb surg30d 
		 flag_apache_above_med  /clodds=pl;
 ods output fitstatistics=fit_red modelanova=model_red cloddspl=clodds_red;

		run;

		/*drop foregin body*/

		
proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race surg30d 
		 flag_apache_above_med  /clodds=pl;
 ods output fitstatistics=fit_red modelanova=model_red cloddspl=clodds_red;

		run;

		/*drop surgery final reduced model*/

proc logistic data=model;
class dialysis_race foreignb surg30d persistent neoplasm 
		steroid flag_apache_above_med flag_age_above_med;
model flag_rsab=dialysis_race 
		 flag_apache_above_med  /clodds=pl;
 ods output fitstatistics=fit_red modelanova=model_red cloddspl=clodds_red;

		run;

proc export data=model_red
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\logistic\model_red.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc export data=clodds_red
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\logistic\clodds_red.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
proc npar1way data=choi_756;
class persistent;
var age;
run;


proc npar1way data=choi_756;
class route;
var age;
run;

proc freq data=model;
tables flag_alive*sab_status/exact;
run;
