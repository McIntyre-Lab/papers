ods listing; ods html;

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/concannon/useful_human_data/aceview_hg19/fusions/sas_data';


/*

I have a set of genes that are expressed

I want to first make a Venn diagram that shows what genes are expressed in what cell types
I already have these

For the set of genes that are expressed in only one cell types I want to flag these
and analyze no further

For genes that are expressed in 2 or more cell types:
For each pairwise comparison (e.g. CD4 vs CD8):
 -- flag genes that have at least 1 exon or one splicing event DE (flag_de)
 -- flag genes that have at least 1 exon or one splicing event DD (flag_as)

Then for each pair, and for genes on in all, count flags


I want to end up with a table:

Cell type expressed (CD4/CD8/CD19, CD4/CD8, etc.)
Genes DD
Genes DS
Genes DE
Genes DS and DE
Genes neither
*/

/* Count the number of DE exons per gene */

data de_exons;
  set con.results_by_fusion_w_flags;
  keep fusion_id flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
run;

data fus2gene;
  set fus.hg19_aceview_fusions_si_info;
  keep fusion_id gene_id;
run;

proc sort data=de_exons;
  by fusion_id;
proc sort data=fus2gene nodup;
  by fusion_id gene_id;
run;

data de_exons2gene;
  merge de_exons (in=in1) fus2gene (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data gene_flags;
   set con.flag_gene_detection_by_cell;
   keep gene_id flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on;
run;

proc sort data=gene_flags;
  by gene_id;
proc sort data=de_exons2gene;
  by gene_id;
run;

data de_exons2gene_2;
  merge gene_flags (in=in1) de_exons2gene (in=in2);
  by gene_id;
  if flag_cd4_gene_on=0 then do;
     flag_cd4cd8_fdr05=.;
     flag_cd4cd19_fdr05=.;
     end;
  if flag_cd8_gene_on=0 then do;
     flag_cd4cd8_fdr05=.;
     flag_cd8cd19_fdr05=.;
     end;
  if flag_cd19_gene_on=0 then do;
     flag_cd4cd19_fdr05=.;
     flag_cd8cd19_fdr05=.;
     end;
run;

proc freq data=de_exons2gene_2 noprint;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05 / out=check;
proc print data=check;
run;

*okay, everything that was tested for DE comes from genes expressed in all three cell types;
* this is good -- I can just count these and merge with the alt splicing flags;

proc sort data=de_exons2gene_2;
  by gene_id;
proc means data=de_exons2gene_2 noprint;
  by gene_id;
  var flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
  output out=exons_de
         sum(flag_cd4cd8_fdr05)=num_exons_de_cd4cd8
         sum(flag_cd4cd19_fdr05)=num_exons_de_cd4cd19
         sum(flag_cd8cd19_fdr05)=num_exons_de_cd8cd19;
run;


/* Count the number of DE events per gene */

data de_events;
   set splicing.flag_splicing_by_gene_dtct_v2;
run;

data de_events_pair;
   set splicing.splicing_results_w_annot_fdr;
   keep event_id flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
run;

proc sort data=de_events;
   by event_id;
proc sort data=de_events_pair;
   by event_id;
run;

data de_events2;
  merge de_events (in=in1) de_events_pair;
  by event_id;
  if in1;
run;

proc freq data=de_events2 noprint;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05 / out=check;
proc print data=check;
run;


/* Need to reset some flags. Drop DE flags if off in that cell type */

data de_events3;
   set de_events2;
   if sum(flag_cd4_gene_on,flag_cd8_gene_on,flag_Cd19_gene_on) < 2 then do;
      flag_cd4cd8_fdr05=.; flag_cd4cd19_fdr05=.; flag_cd8cd19_fdr05=.; flag_anova_fdr_05=.;
      if flag_Cd4_gene_on=0 then flag_cd4_on=0;
      if flag_Cd8_gene_on=0 then flag_cd8_on=0;
      if flag_Cd19_gene_on=0 then flag_cd19_on=0;
      end;
   else do;
     if flag_cd4_gene_on=0 then do;
         flag_cd4cd8_fdr05=.; flag_cd4cd19_fdr05=.;
         flag_cd4_on=0;
         end;
     if flag_cd8_gene_on=0 then do;
         flag_cd4cd8_fdr05=.; flag_cd8cd19_fdr05=.;
         flag_cd8_on=0;
         end;
     if flag_cd19_gene_on=0 then do;
         flag_cd4cd19_fdr05=.; flag_cd8cd19_fdr05=.;
         flag_cd19_on=0;
         end;
     end;
run;

proc sort data=de_events3;
  by gene_id;
proc means data=de_events3 noprint;
  by gene_id;
  var  flag_cd4_on flag_cd8_on flag_cd19_on
       flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05 flag_anova_fdr_05;
  output out=events_de
         sum(flag_cd4_on)=num_events_on_cd4
         sum(flag_cd8_on)=num_events_on_cd8
         sum(flag_cd19_on)=num_events_on_cd19
         sum(flag_anova_fdr_05)=num_events_de_overall
         sum(flag_cd4cd8_fdr05)=num_events_de_cd4cd8
         sum(flag_cd4cd19_fdr05)=num_events_de_cd4cd19
         sum(flag_cd8cd19_fdr05)=num_events_de_cd8cd19;
run;

proc sort data=events_de;
  by gene_id;
proc sort data=exons_de;
  by gene_id;
run;

data genes_de;
  merge exons_de events_de;
  by gene_id;
  drop _TYPE_ _FREQ_;
run;


* Summary of alternative splicing -- this has the DD and DS pieces;
data diff_splicing_summary;
   set con.gene_diff_splicing_summary;
run;

proc sort data=genes_de;
  by gene_id;
proc sort data=diff_splicing_summary;
  by gene_id;
run;

data diff_splicing;
  merge diff_splicing_summary (in=in1) genes_de (in=in2);
  by gene_id;
  if in1;
run;

/* For each gene, I want to flag if it is DD, DS and DE for each pairwise comparison */

data flag_dd_ds_de;
   set diff_splicing;
   /* DD flag is kind of already set: this is the "gene_on" flags and "flag_gene_on_multicell" */

   if flag_gene_on_multicell=1 then do;

     /* Fix some counts */
     if flag_cd4_gene_on=0 then do;
        num_exons_de_cd4cd8=0; num_events_de_cd4cd8=0; 
        num_exons_de_cd4cd19=0; num_events_de_cd4cd19=0; end;

     if flag_cd8_gene_on=0 then do;
        num_exons_de_cd4cd8=0; num_events_de_cd4cd8=0; 
        num_exons_de_cd8cd19=0; num_events_de_cd8cd19=0; end;

     if flag_cd19_gene_on=0 then do;
        num_exons_de_cd4cd19=0; num_events_de_cd4cd19=0; 
        num_exons_de_cd8cd19=0; num_events_de_cd8cd19=0; end;


     /* Flag if gene is DS for each comparison */
     if flag_cd4cd8_exon_dd=1 or flag_cd4cd8_event_dd=1 then flag_cd4cd8_gene_ds=1;
     else flag_cd4cd8_gene_ds=0;
     if flag_cd4cd19_exon_dd=1 or flag_cd4cd19_event_dd=1 then flag_cd4cd19_gene_ds=1;
     else flag_cd4cd19_gene_ds=0;
     if flag_cd8cd19_exon_dd=1 or flag_cd8cd19_event_dd=1 then flag_cd8cd19_gene_ds=1;
     else flag_cd8cd19_gene_ds=0;

     /* Flag if gene is DE for each comparison */
     if num_exons_de_cd4cd8 > 0 or num_events_de_cd4cd8 > 0 then flag_cd4cd8_gene_de=1;
     else flag_cd4cd8_gene_de=0;
     if num_exons_de_cd4cd19 > 0 or num_events_de_cd4cd19 > 0 then flag_cd4cd19_gene_de=1;
     else flag_cd4cd19_gene_de=0;
     if num_exons_de_cd8cd19 > 0 or num_events_de_cd8cd19 > 0 then flag_cd8cd19_gene_de=1;
     else flag_cd8cd19_gene_de=0;

     if flag_cd4cd8_gene_ds=1 or flag_cd4cd19_gene_ds=1 or flag_cd8cd19_gene_ds=1 then flag_gene_ds=1;
     else flag_gene_ds=0;

     if flag_cd4cd8_gene_de=1 or flag_cd4cd19_gene_de=1 or flag_cd8cd19_gene_de=1 then flag_gene_de=1;
     else flag_gene_de=0;

   end;
run;


/* Make permenant so I can refer to this later */

data con.flag_genes_dd_ds_de;
  set flag_dd_ds_de;
run;

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=flag_dd_ds_de nodup;
   by gene_id;
proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d nodup;
   by gene_id;
run;


data flag_dd_ds_de_immuno;
  merge flag_dd_ds_de (in=in1) ai (in=in2) t1d (in=in3);
  by gene_id;
  if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
  if in1 then output;
run;
   

/* COUNT!!! */

proc freq data=flag_dd_ds_de_immuno noprint;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*
          flag_cd4cd8_gene_ds*flag_cd4cd19_gene_ds*flag_cd8cd19_gene_ds*
          flag_cd4cd8_gene_de*flag_cd4cd19_gene_de*flag_cd8cd19_gene_de / out = gene_counts;
run;

proc freq data=flag_dd_ds_de_immuno noprint;
   where flag_immuno_gene=1;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*
          flag_cd4cd8_gene_ds*flag_cd4cd19_gene_ds*flag_cd8cd19_gene_ds*
          flag_cd4cd8_gene_de*flag_cd4cd19_gene_de*flag_cd8cd19_gene_de / out = ai_gene_counts;
run;

proc freq data=flag_dd_ds_de_immuno noprint;
   where flag_immunobase_diabetes_gene=1;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*
          flag_cd4cd8_gene_ds*flag_cd4cd19_gene_ds*flag_cd8cd19_gene_ds*
          flag_cd4cd8_gene_de*flag_cd4cd19_gene_de*flag_cd8cd19_gene_de / out = t1d_gene_counts;
run;


proc export data=gene_counts
     outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/dd_de_ds_gene_counts2_all2.csv"
     dbms=csv replace;
run;
proc export data=ai_gene_counts
     outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/dd_de_ds_gene_counts2_autoimmune2.csv"
     dbms=csv replace;
run;
proc export data=t1d_gene_counts
     outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/dd_de_ds_gene_counts2_diabetes2.csv"
     dbms=csv replace;
run;



proc freq data=flag_dd_ds_de_immuno noprint;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*
          flag_gene_ds*flag_gene_de / out = gene_counts;
run;

proc freq data=flag_dd_ds_de_immuno noprint;
   where flag_immuno_gene=1;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*
          flag_gene_ds*flag_gene_de / out = ai_gene_counts;
run;

proc freq data=flag_dd_ds_de_immuno noprint;
   where flag_immunobase_diabetes_gene=1;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*
          flag_gene_ds*flag_gene_de / out = t1d_gene_counts;
run;

proc export data=gene_counts
     outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/dd_de_ds_gene_counts_simple_all2.csv"
     dbms=csv replace;
run;

proc export data=ai_gene_counts
     outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/dd_de_ds_gene_counts_simple_autoimmune2.csv"
     dbms=csv replace;
run;

proc export data=t1d_gene_counts
     outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/dd_de_ds_gene_counts_simple_diabetes2.csv"
     dbms=csv replace;
run;

