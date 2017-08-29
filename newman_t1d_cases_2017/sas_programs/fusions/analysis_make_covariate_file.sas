
/* Make covariate file for running models */

libname con '/home/jrbnewman/concannon/sas_data';

data peer_factors;
  set con.peer_factors;
run;

data covariate_data;
  set con.design_by_subject_new;
run;

proc sort data=peer_factors;
  by Name;
proc sort data=covariate_data;
  by Name;
run;

data all_covariates;
   merge covariate_data (in=in1) peer_factors (in=in2);
   by Name;
   if in1 and in2 then output;
run;

data con.all_covariates;
   set all_covariates;
run;
