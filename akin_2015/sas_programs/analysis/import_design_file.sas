/* Libraries */

libname sugrue '!MCLAB/sugrue/sas_data';

/* import sample key */

proc import datafile='!MCLAB/sugrue/design_files/sample_key.csv'
   out=sugrue.sample_key
   dbms=csv
   replace;
   guessingrows=100;
run;


proc import datafile='!MCLAB/sugrue/design_files/subject2sample.csv'
   out=sugrue.subject2sample
   dbms=csv
   replace;
   guessingrows=100;
run;
