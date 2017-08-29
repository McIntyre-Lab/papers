/* Import autoimmune candidate gene flags */

libname con '/home/jrbnewman/concannon/sas_data';

proc import datafile='/home/jrbnewman/concannon/design_files/immunogene_flags.csv' out=immunogene_flags
dbms=csv replace; guessingrows=73000;
run;

data con.immunogene_flags;
   set immunogene_flags;
run;
