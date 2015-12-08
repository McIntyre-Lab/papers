/* Get splicing event 2 gene_id prepared for HPC */

libname splice '/mnt/data/splice';
libname eqtl '/mnt/data/eqtls/sas_data';

data splicing2gene;
   set splice.splicing_events_annotations;
   keep event_id gene_id;
run;

data eqtl.splicing2gene;
   set splicing2gene;
run;

