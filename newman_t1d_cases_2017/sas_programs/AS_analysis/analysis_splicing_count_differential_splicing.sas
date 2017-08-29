
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


/* Merge gene-level flags with "cleaned" splicing events */

data splicing_flags;
   set splicing.splicing_results_clean;
run;

proc sort data=gene_flags;
   by gene_id;
proc sort data=splicing_flags;
   by gene_id;
run;

data splicing_gene_flags;
  merge splicing_flags (in=in1) gene_flags (in=in2);
  by gene_id;
  if in1 and in2;
run;


/* Flag events as: cell-specific (ie, gene is only detected in one cell type)
                    alternatively spliced (gene is in at least two, but event is in fewer) */

data flag_splicing_spec;
   set splicing_gene_flags;
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
        flag_anova_fdr_05=.;
        end;
   if flag_cd8_gene_on=0 then do;
        flag_cd8_on=0;
        flag_anova_fdr_05=.;
        end;
   if flag_cd19_gene_on=0 then do;
        flag_cd19_on=0;
        flag_anova_fdr_05=.;
        end;
   sum_flags=flag_cd4_on+flag_cd8_on+flag_cd19_on;
run;




/* Counts! */

%macro counts(prefix,outname);
proc freq data=flag_splicing_spec noprint;
   %if &prefix. = ja %then %do; where flag_junction_annotated=1; %end;
   %if &prefix. = ju %then %do; where flag_junction_annotated=0 and flag_intron_retention=0; %end;
   %if &prefix. = ir %then %do; where flag_intron_retention=1; %end;
   %if &prefix. = es %then %do; where flag_exonskip=1; %end;
   %if &prefix. = ad %then %do; where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0; %end;
   %if &prefix. = aa %then %do; where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0; %end;
   %if &prefix. = ada %then %do; where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0; %end;
   tables flag_cell_specific*flag_alt_spliced*flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4_on*flag_cd8_on*flag_cd19_on*flag_anova_fdr_05*sum_flags / out = &prefix._counts;
run;



/* Export data */

data &prefix._counts2;
   set &prefix._counts;
   drop PERCENT;
run;


/* Export data */

proc export data=&prefix._counts2
     outfile="!PATCON/pipeline_output/differential_splicing_counts/counts_by_event_and_gene_detection_&outname..csv"
     dbms=csv replace;
run;


%mend;

%counts(all,all_events);
%counts(ja,annotated_junctions);
%counts(ju,unannotated_junctions);
%counts(ir,intron_retention);
%counts(es,exonskipping);
%counts(ad,alt_donor);
%counts(aa,alt_acceptor);
%counts(ada,alt_donor_and_acceptor);

/* Make permenant */

data splicing.flag_splicing_by_gene_dtct_v2;
   set flag_splicing_spec;
run;


