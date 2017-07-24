ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For detected unannotated splicing events (unannotated junctions, IR events), infer what genes may have expressed novel isoforms.

Counting the number of detected unannotated junctions and IR events per gene. This gives an idea as to how many genes have potentially novel transcripts (albeit, we cannot determine what they are from just short read data).

I am going to limit this to ONLY expressed genes */

data exp_genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

data detected_events;
   set event.splicing_on_apn_gt0;
   where flag_junction_annotated=0 and flag_splicing_on=1;
   keep event_id flag_junction_annotated flag_intron_retention;
run;

data event2gene;
   set evspl.splicing_Events_annot_refseq;
   keep event_id gene_id;
run;

proc sort data=detected_events;
   by event_id;
proc sort data=event2gene;
   by event_id;
run;

data dtct_event2gene;
   merge event2gene (in=in1) detected_events (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=dtct_event2gene;
   by gene_id;
proc sort data=exp_genes;
   by gene_id;
run;

data dtct_novel_exp_genes;
   merge exp_genes (in=in1) dtct_event2gene (in=in2);
   by gene_id;
   if in1 and in2;
run;

/* Count unannotated events per gene */

* Unannotated junctions;
proc freq data=dtct_novel_exp_genes noprint;
  where flag_intron_retention=0;
  tables gene_id / out=unannot_by_gene;
run;

* IR events;
proc freq data=dtct_novel_exp_genes noprint;
  where flag_intron_retention=1;
  tables gene_id / out=ir_by_gene;
run;

data unannot_by_gene2;
   set unannot_by_gene;
   rename count=num_unannotated_junctions;
   drop percent;
run;

data ir_by_gene2;
   set ir_by_gene;
   rename count=num_ir_events;
   drop percent;
run;

proc sort data=unannot_by_gene2;
  by gene_id;
proc sort data=ir_by_gene2;
  by gene_id;
run;

data novel_by_gene;
  merge unannot_by_gene2 (in=in1) ir_by_gene2 (in=in2);
  by gene_id;
  if not in1 then num_unannotated_junctions=0;
  if not in2 then num_ir_events=0;
  if num_unannotated_junctions > 0 then flag_gene_has_unannotated_junc=1;
  else  flag_gene_has_unannotated_junc=0;
  if  num_ir_events > 0 then flag_gene_has_ir=1;
  else  flag_gene_has_ir=0;
run;

proc freq data=novel_by_gene;
  tables flag_gene_has_unannotated_junc*flag_gene_has_ir;
run;


/*
   flag_gene_has_unannotated_junc
             flag_gene_has_ir

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |   7425 |   7425
            |   0.00 |  55.90 |  55.90
            |   0.00 | 100.00 |
            |   0.00 |  58.61 |
   ---------+--------+--------+
          1 |    613 |   5244 |   5857
            |   4.62 |  39.48 |  44.10
            |  10.47 |  89.53 |
            | 100.00 |  41.39 |
   ---------+--------+--------+
   Total         613    12669    13282
                4.62    95.38   100.00

7425 genes with putative intron retention
613 genes with an unannotated junction
5244 genes with putative IR and unannotated junctions
*/

/* Make permenant */

data event.unannotated_events_by_gene;
   set novel_by_gene;
run;

	
