/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Making a big counts table. Per gene:
   flag if cell type specific (expressed in only one cell type)
   for genes in at least two cell types:
	- exons DD in CD4vsCD8, CD4vsCD19, CD8vsCD19
	- events DD in CD4vsCD8, CD4vsCD19, CD8vsCD19
	- flag if any different!

   Then calculate the %genes with any different. This will be the % genes alternatively spliced */

/* Get gene detection flags */

data gene_flags; *only keep genes expressed in at least one cell type;
  set con.flag_gene_detection_by_cell;
  if flag_cd4_gene_on+flag_cd8_gene_on+flag_cd19_gene_on > 1
       then flag_gene_on_multicell=1;
       else flag_gene_on_multicell=0;
  if flag_cd4_gene_on+flag_cd8_gene_on+flag_cd19_gene_on=0 then delete;
  keep gene_id flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on flag_gene_on_multicell;
run; *47063 genes;

/* Get fusion detection flags */

data fusion_flags;
   set con.fusions_on_gt_apn0;
   /* CD4 vs CD8 differential detection */
   if flag_cd4_on=1 and flag_cd8_on=0 then flag_cd4cd8_exon_dd=1;
   else if flag_cd4_on=0 and flag_cd8_on=1 then flag_cd4cd8_exon_dd=1;
   else flag_cd4cd8_exon_dd=0;
   /* CD4 vs CD19 differential detection */
   if flag_cd4_on=1 and flag_cd19_on=0 then flag_cd4cd19_exon_dd=1;
   else if flag_cd4_on=0 and flag_cd19_on=1 then flag_cd4cd19_exon_dd=1;
   else flag_cd4cd19_exon_dd=0;
   /* CD8 vs CD19 differential detection */
   if flag_cd8_on=1 and flag_cd19_on=0 then flag_cd8cd19_exon_dd=1;
   else if flag_cd8_on=0 and flag_cd19_on=1 then flag_cd8cd19_exon_dd=1;
   else flag_cd8cd19_exon_dd=0;
   keep fusion_id flag_cd4cd8_exon_dd flag_cd4cd19_exon_dd flag_cd8cd19_exon_dd;
run;

/* Get fusion2gene */
data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

proc sort data=fus2gene;
  by fusion_id;
proc sort data=fusion_flags;
  by fusion_id;
run;

data fusion_flags_w_gene;
   merge fus2gene (in=in1) fusion_flags (in=in2);
   by fusion_id;
   if in1 and in2;
run;

proc sort data=fusion_flags_w_gene;
   by gene_id;
proc means data=fusion_flags_w_gene noprint;
   by gene_id;
   var flag_cd4cd8_exon_dd flag_cd4cd19_exon_dd flag_cd8cd19_exon_dd;
   output out=num_fus_dd_by_gene
          sum(flag_cd4cd8_exon_dd)=num_exons_dd_cd4cd8
          sum(flag_cd4cd19_exon_dd)=num_exons_dd_cd4cd19
          sum(flag_cd8cd19_exon_dd)=num_exons_dd_cd8cd19;
run;


data flag_fus_dd_by_gene;
   set num_fus_dd_by_gene;
   if num_exons_dd_cd4cd8 > 0 then flag_cd4cd8_exon_dd=1; else flag_cd4cd8_exon_dd=0;
   if num_exons_dd_cd4cd19 > 0 then flag_cd4cd19_exon_dd=1; else flag_cd4cd19_exon_dd=0;
   if num_exons_dd_cd8cd19 > 0 then flag_cd8cd19_exon_dd=1; else flag_cd8cd19_exon_dd=0;
run;

/* Get event detection flags */

data event_flags;
   set splicing.splicing_results_clean;
   /* CD4 vs CD8 differential detection */
   if flag_cd4_on=1 and flag_cd8_on=0 then flag_cd4cd8_event_dd=1;
   else if flag_cd4_on=0 and flag_cd8_on=1 then flag_cd4cd8_event_dd=1;
   else flag_cd4cd8_event_dd=0;
   /* CD4 vs CD19 differential detection */
   if flag_cd4_on=1 and flag_cd19_on=0 then flag_cd4cd19_event_dd=1;
   else if flag_cd4_on=0 and flag_cd19_on=1 then flag_cd4cd19_event_dd=1;
   else flag_cd4cd19_event_dd=0;
   /* CD8 vs CD19 differential detection */
   if flag_cd8_on=1 and flag_cd19_on=0 then flag_cd8cd19_event_dd=1;
   else if flag_cd8_on=0 and flag_cd19_on=1 then flag_cd8cd19_event_dd=1;
   else flag_cd8cd19_event_dd=0;
   keep event_id flag_cd4cd8_event_dd flag_cd4cd19_event_dd flag_cd8cd19_event_dd;
run;

/* Get splicing2gene */
data event2gene;
   set splice.splicing_events_annotations;
   keep event_id gene_id;
run;

proc sort data=event2gene;
  by event_id;
proc sort data=event_flags;
  by event_id;
run;

data event_flags_w_gene;
   merge event2gene (in=in1) event_flags (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=event_flags_w_gene;
   by gene_id;
proc means data=event_flags_w_gene noprint;
   by gene_id;
   var flag_cd4cd8_event_dd flag_cd4cd19_event_dd flag_cd8cd19_event_dd;
   output out=num_event_dd_by_gene
          sum(flag_cd4cd8_event_dd)=num_events_dd_cd4cd8
          sum(flag_cd4cd19_event_dd)=num_events_dd_cd4cd19
          sum(flag_cd8cd19_event_dd)=num_events_dd_cd8cd19;
run;


data flag_event_dd_by_gene;
   set num_event_dd_by_gene;
   if num_events_dd_cd4cd8 > 0 then flag_cd4cd8_event_dd=1; else flag_cd4cd8_event_dd=0;
   if num_events_dd_cd4cd19 > 0 then flag_cd4cd19_event_dd=1; else flag_cd4cd19_event_dd=0;
   if num_events_dd_cd8cd19 > 0 then flag_cd8cd19_event_dd=1; else flag_cd8cd19_event_dd=0;
run;

/* Merge with gene-level flags */

proc sort data=flag_event_dd_by_gene;
   by gene_id;
proc sort data=flag_fus_dd_by_gene;
   by gene_id;
proc sort data=gene_flags;
   by gene_id;
run;

data gene_dd_summary;
   merge gene_flags (in=in1) flag_fus_dd_by_gene (in=in2)
         flag_event_dd_by_gene (in=in3);
   by gene_id;
   if not in2 then do;
   flag_cd4cd8_exon_dd=0;
   flag_cd4cd19_exon_dd=0;
   flag_cd8cd19_exon_dd=0;
   num_exons_dd_cd4cd8=0;
   num_exons_dd_cd4cd19=0;
   num_exons_dd_cd8cd19=0;
   end;
   if not in3 then do;
   flag_cd4cd8_event_dd=0;
   flag_cd4cd19_event_dd=0;
   flag_cd8cd19_event_dd=0;
   num_exons_dd_cd4cd8=0;
   num_exons_dd_cd4cd19=0;
   num_exons_dd_cd8cd19=0;
   end;
   if in1 then output;
run;

data gene_dd_summary2;
  set gene_dd_summary;
  if flag_gene_on_multicell=0 then do;
     *set these to missing as these are genes only on in one cell type;
     flag_cd4cd8_exon_dd=.;
     flag_cd4cd19_exon_dd=.;
     flag_cd8cd19_exon_dd=.;
     flag_cd4cd8_event_dd=.;
     flag_cd4cd19_event_dd=.;
     flag_cd8cd19_event_dd=.;
     num_exons_dd_cd4cd8=.;
     num_exons_dd_cd4cd19=.;
     num_exons_dd_cd8cd19=.;
     num_events_dd_cd4cd8=.;
     num_events_dd_cd4cd19=.;
     num_events_dd_cd8cd19=.;
     end;
  else do;
     if flag_cd4_gene_on=0 then do; *i.e. if gene is off, then set comparisons to 0;
          num_exons_dd_cd4cd8=.; num_exons_dd_cd4cd19=.;
          flag_cd4cd8_exon_dd=0; flag_cd4cd19_exon_dd=0;
          num_events_dd_cd4cd8=.; num_events_dd_cd4cd19=.;
          flag_cd4cd8_event_dd=0; flag_cd4cd19_event_dd=0;
          end;
     if flag_cd8_gene_on=0 then do; *i.e. if gene is off, then set comparisons to 0;
          num_exons_dd_cd4cd8=.; num_exons_dd_cd8cd19=.;
          flag_cd4cd8_exon_dd=0; flag_cd8cd19_exon_dd=0;
          num_events_dd_cd4cd8=.; num_events_dd_cd8cd19=.;
          flag_cd4cd8_event_dd=0; flag_cd8cd19_event_dd=0;
          end;
     if flag_cd19_gene_on=0 then do; *i.e. if gene is off, then set comparisons to 0;
          num_exons_dd_cd4cd19=.; num_exons_dd_cd8cd19=.;
          flag_cd4cd19_exon_dd=0; flag_cd8cd19_exon_dd=0;
          num_events_dd_cd4cd19=.; num_events_dd_cd8cd19=.;
          flag_cd4cd19_event_dd=0; flag_cd8cd19_event_dd=0;
          end;
     if flag_cd4cd8_exon_dd=1 or flag_cd4cd8_event_dd=1 or 
        flag_cd4cd19_exon_dd=1 or flag_cd4cd19_event_dd=1 or 
        flag_cd8cd19_exon_dd=1 or flag_cd8cd19_event_dd=1 
     then flag_gene_any_dd=1;
     else flag_gene_any_dd=0;
     end;
run;

proc freq data=gene_dd_summary2;
   tables flag_gene_any_dd;
run;

/*

                                              Cumulative    Cumulative
 flag_gene_any_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       28961       65.84         28961        65.84
                1       15024       34.16         43985       100.00

*/

/* Make permenant */

data con.gene_diff_splicing_summary;
  set gene_dd_summary2;
run;


