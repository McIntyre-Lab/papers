libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";

/* Michelle would like the final results file with the enrichment test results
 * combined. */

proc sort data=ribo.motif_enrichment_revision;
   by primary_fbgn;
   run;
proc sort data=ribo.motif_flags_and_cnts;
  by primary_fbgn;
  run;

data motif_enrichment_flags_revision;
  merge ribo.motif_enrichment_revision (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1;
  run;


proc export data=motif_enrichment_flags_revision
    outfile="!MCLAB/arbeitman/arbeitman_ribotag/reports/gene_list_with_enrichment_revision.csv"
    DBMS=CSV REPLACE;
    run;
