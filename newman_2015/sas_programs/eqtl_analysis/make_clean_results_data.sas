/* Make a clean eQTL summary dataset */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';

/* Get list of cleaned splicing events */

data splicing_clean;
   set splicing.splicing_results_clean;
   keep event_id;
   rename event_id=feature_id;
run;

/* Split eQTL results */

data eqtl_exons eqtl_splicing;
   set eqtl.eqtl_results_summary;
   if feature_type='exon' then output eqtl_exons;
   else output eqtl_splicing;
run;

proc sort data=splicing_clean;
   by feature_id;
proc sort data=eqtl_splicing;
   by feature_id;
run;

data eqtl_splicing_clean no_exp no_eqtl;
   merge eqtl_splicing (in=in1) splicing_clean (in=in2);
   by feature_id;
   if in1 and in2 then output eqtl_splicing_clean;
   else if in1 then outpot no_exp;
   else output no_eqtl;
run;

data eqtl_results_clean;
   set eqtl_exons eqtl_splicing_clean;
run;

/* Make permenant */

data eqtl.eqtl_results_summary_clean;
   set eqtl_results_clean;
run;
