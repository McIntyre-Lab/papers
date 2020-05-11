
/* data from Recur SAB_isolates_Data_Codebook_2018_05_18_Lauren3-17-20*/
/*sheet IsolatedCheckedOut*/

proc import out=work.pfge
 datafile="Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\relapse_data_from_duke\pfge_4_sas.csv"
 dbms=csv replace;
 getnames=yes;
 datarow=2;
 run;

proc freq data=pfge noprint;
tables patientnumber/out=count_patient;
run;

proc freq data=pfge;
tables PFGE;
run;

/* fix coding for the one patient with the PFGE number still included
3834	65	23-24 ID
3850	65	23-24 ID
3899	65	25-26 ID
4122	65	25-26 ID */

data pfge_fixed;
set pfge;
/*recode below set for patient 65 vance verified email*/
if id=3834 then pfge="Identical";
if id=3850 then pfge="Different";
if id=3899 then pfge="Identical";
if id=4122 then pfge="Identical";
 /*fix typos*/
if PFGE="identical" then PFGE="Identical";
if PFGE="different" then PFGE="Different";

/*fix days for patient 35 intake was 5/3/11 and second episode was 7/8/2011*/

if id=4697 then days=66;
if id=4739 then days="";

/*the last episode has a code, this does not make sense in the format 
where the previous epsiode has the sequential value for the next episode so remove to avoid later confusion*/
if days= . then PFGE="";

/* patient 83 has 4 episodes but one was excluded called identical.*/
/*felicia verified episode 3838 was a double count on intake and without sample this numbers match*/
if id=3838 then delete;


/*patient 0 is complicated see Felicia's email, to avoid recoding everything I am putting days with the later case*/

if id=4130 then days=371;
if id=4130 then pfge="Different";
if id=4132 then days=.;
if id=4132 then pfge="";

run;




/* this is 72 remove excluded patients*/
data choi.pfge_included;
set pfge_fixed;
if patientnumber=63 or patientnumber=68 or patientnumber=71 then delete;
run;

/*154 episodes 69 patients*/

proc freq data=pfge_included;
tables PFGE;
run;

proc freq data=pfge_included noprint;
tables patientnumber/out=count_patient;
run;
/*this has 69*/
