
/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Determine the set of genes that are detected. Then for the set of genes detected in all three cell types,
   count the number "on". The number that are DD are now evidence of alternative splicing!

   I would like to count the following:

   cd4*cd8*cd19 for events from genes on in all cell types

   cd4*cd8 for events from genes only on in CD4 and CD8
   cd4*cd19 for events from genes only on in CD4 and CD19
   cd8*cd19 for events from genes only on in CD8 and CD19

   cd4 events from genes only on in CD4
   cd8 events from genes only on in CD8
   cd19 events from genes only on in CD19

*/

/* Get gene-level flags */

data gene_flags;
   set con.flag_gene_detection_by_cell;
   keep gene_id flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on;
run;


/* Get gene to fusion: will need to collapse these */

data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

proc sort data=fus2gene;
  by gene_id;
proc sort data=gene_flags;
  by gene_id;
run;

data fus2gene_flags;
  merge fus2gene (in=in1) gene_flags (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=fus2gene_flags;
  by fusion_id;
proc means data=fus2gene_flags noprint;
  by fusion_id;
  var flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on;
  output out=fus_w_gene_flags max=;
run;



/* Fusion on flags */
data on_flags;
  set con.fusions_on_gt_apn0;
  drop flag_fusion_on0;
run;

data results;
  set con.results_by_fusion_final;
  keep fusion_id flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
run;


/* Merge gene-level flags with fusion data */


proc sort data=fus_w_gene_flags;
  by fusion_id;
proc sort data=on_flags;
  by fusion_id;
proc sort data=results;
  by fusion_id;
run;


data fusion_gene_flags;
  merge fus_w_gene_flags (in=in1) on_flags (in=in2) results (in=in3);
  by fusion_id;
run;

/* Flag fusions as: cell-specific (ie, gene is only detected in one cell type)
                    alternatively spliced (gene is in at least two, but fusion is in fewer) */

data flag_fusion_spec;
   set fusion_gene_flags;
   /* Delete fusions that are off in all cell-types */
   if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=0 then delete;

   /* Flag cell-specific fusions */
   if flag_cd4_gene_on=1 and flag_cd8_gene_on=0 and flag_cd19_gene_on=0 then flag_cell_specific=1;
   else if flag_cd4_gene_on=0 and flag_cd8_gene_on=1 and flag_cd19_gene_on=0 then flag_cell_specific=1;
   else if flag_cd4_gene_on=0 and flag_cd8_gene_on=0 and flag_cd19_gene_on=1 then flag_cell_specific=1;
   else flag_cell_specific=0;

   /* Flag alternatively spliced fusions */
   if flag_cd4_gene_on=1 and flag_cd8_gene_on=1 and flag_cd19_gene_on=0 then do;
      if flag_cd4_on=1 and flag_cd8_on=0 then flag_alt_spliced=1;
      else if flag_cd4_on=0 and flag_cd8_on=1 then flag_alt_spliced=1;
      else flag_alt_spliced=0;
   end;
   else if flag_cd4_gene_on=1 and flag_cd8_gene_on=0 and flag_cd19_gene_on=1 then do;
      if flag_cd4_on=1 and flag_cd19_on=0 then flag_alt_spliced=1;
      else if flag_cd4_on=0 and flag_cd19_on=1 then flag_alt_spliced=1;
      else flag_alt_spliced=0;
   end;
   else if flag_cd4_gene_on=0 and flag_cd8_gene_on=1 and flag_cd19_gene_on=1 then do;
      if flag_cd8_on=1 and flag_cd19_on=0 then flag_alt_spliced=1;
      else if flag_cd8_on=0 and flag_cd19_on=1 then flag_alt_spliced=1;
      else flag_alt_spliced=0;
   end;
   else if flag_cd4_gene_on=1 and flag_cd8_gene_on=1 and flag_cd19_gene_on=1 then do;
      if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) < 3 then flag_alt_spliced=1;
      else flag_alt_spliced=0;
   end;
   else flag_alt_spliced=0;

   /* Reset fusion on flags if gene is off */
   if flag_cd4_gene_on=0 then do;
        flag_cd4_on=0;
        flag_cd4cd8_fdr05=.;
        flag_cd4cd19_fdr05=.;
        end;
   if flag_cd8_gene_on=0 then do;
        flag_cd8_on=0;
        flag_cd4cd8_fdr05=.;
        flag_cd8cd19_fdr05=.;
        end;
   if flag_cd19_gene_on=0 then do;
        flag_cd19_on=0;
        flag_cd8cd19_fdr05=.;
        flag_cd4cd19_fdr05=.;
        end;
   sum_flags=flag_cd4_on+flag_cd8_on+flag_cd19_on;
run;



/* Counts! */

proc freq data=flag_fusion_spec noprint;
   tables flag_cell_specific*flag_alt_spliced*
          flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4_on*flag_cd8_on*flag_cd19_on*
          flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05*sum_flags / out = fusion_counts;
run;

/* Make dataset permenant */

data con.flag_fusions_by_gene_detection2;
    set flag_fusion_spec;
run;


/* Export data */

data fusion_counts2;
   set fusion_counts;
   drop PERCENT;
run;


proc export data=fusion_counts2
     outfile="!PATCON/pipeline_output/differential_splicing_counts/counts_by_fusion_and_gene_detection.csv"
     dbms=csv replace;
run;

