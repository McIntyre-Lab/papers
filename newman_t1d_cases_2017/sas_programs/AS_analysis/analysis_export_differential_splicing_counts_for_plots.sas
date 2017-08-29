/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* I want to count the number of DD events per gene and plot */

data diff_splicing;
   set con.gene_diff_splicing_summary;
   /* Count total number of DD exons */

   num_exons_dd_total=sum(num_exons_dd_cd4cd8,num_exons_dd_cd4cd19,num_exons_dd_cd8cd19);

   /* Count total number of DD events */
   num_events_dd_total=sum(num_events_dd_cd4cd8,num_events_dd_cd4cd19,num_events_dd_cd8cd19);

   /* Count total number of DD features */
   num_features_dd_total=sum(num_exons_dd_total,num_events_dd_total);
   if num_features_dd_total=. and flag_gene_on_multicell=1 then num_features_dd_total=0;
   drop _TYPE_ _FREQ_;
run;

/* Export data for plots */

data data_for_export;
   set diff_splicing;
   where flag_gene_on_multicell=1;
run;

proc export data=data_for_export
     outfile="!PATCON/pipeline_output/differential_splicing_counts/differential_splicing_counts_for_plots.csv"
     dbms=csv replace;
run;

