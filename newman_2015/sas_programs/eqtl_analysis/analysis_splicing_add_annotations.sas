/* Merge in splicing annotations  (event2gene) */

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data2'; 
libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';

/* Merge splicing counts into one file */

data mysas.splicing_counts_for_eqtls;
set mysas.splicing_counts_for_eqtls_1 mysas.splicing_counts_for_eqtls_2 mysas.splicing_counts_for_eqtls_3 mysas.splicing_counts_for_eqtls_4
    mysas.splicing_counts_for_eqtls_5 mysas.splicing_counts_for_eqtls_6 mysas.splicing_counts_for_eqtls_7 mysas.splicing_counts_for_eqtls_8
    mysas.splicing_counts_for_eqtls_9 mysas.splicing_counts_for_eqtls_10 mysas.splicing_counts_for_eqtls_11 mysas.splicing_counts_for_eqtls_12
    mysas.splicing_counts_for_eqtls_13 mysas.splicing_counts_for_eqtls_14 mysas.splicing_counts_for_eqtls_15 mysas.splicing_counts_for_eqtls_16
    mysas.splicing_counts_for_eqtls_17 mysas.splicing_counts_for_eqtls_18 mysas.splicing_counts_for_eqtls_19 mysas.splicing_counts_for_eqtls_20
    mysas.splicing_counts_for_eqtls_21 mysas.splicing_counts_for_eqtls_22 mysas.splicing_counts_for_eqtls_23 mysas.splicing_counts_for_eqtls_24
    mysas.splicing_counts_for_eqtls_25 ;
run;


/* Merge gene ids */

proc sort data=eqtl.splicing2gene;
   by event_id;
proc sort data=mysas.splicing_counts_for_eqtls;
   by event_id;
run;

data eqtl.splicing_counts no_gene;
   merge eqtl.splicing2gene (in=in1) mysas.splicing_counts_for_eqtls (in=in2);
   by event_id;
   if in1 and in2 then output eqtl.splicing_counts;
   else if in2 then output no_gene;
run;

