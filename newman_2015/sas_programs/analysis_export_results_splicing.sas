libname splicing '/mnt/data/splicing/';

/* Export splicing results to Concannon share */


proc export data=splicing.splicing_results_w_annot outfile="/home/jrbnewman/concannon/reports/splicing/results_by_splicing.csv"
   dbms=csv replace;
run;


