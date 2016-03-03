/* Import design file into SAS */

/* Libraries */
libname con '/home/jrbnewman/concannon/sas_data'; *local/share libname;

proc import datafile='/home/jrbnewman/concannon/design_files/diabetes_design_file_pe_covar.csv'
   out=design_by_sample dbms=tab replace; guessingrows=5024;
run;

/* Make permenant */

data con.design_by_subject_new;
   set design_by_sample;
run;


