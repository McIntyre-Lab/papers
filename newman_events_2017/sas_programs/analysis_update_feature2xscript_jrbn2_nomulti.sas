ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Update the assignment of features to genes and transcripts by removing transcripts of genes
   that are not expressed, transcripts for genes that only have multigene exonic regions (too tricky to handle)
   and any transcript of expressed genes that has no features detected. For this latter category, I am going to
   put these transcripts into a separate list for exploring further */


* xscripts of expressed genes;
data xs_gene_no_exp;
  set event.flag_xscript_w_gene_on;
  if flag_xscript_gene_exp=0 or flag_gene_only_mult=1;
  keep transcript_id;
run;


proc sort data=xs_gene_no_exp;
  by transcript_id;
run;

/* get list xscripts without detected features */

* Feature2xscript list;
data feat2xs2gene;
  set event.feature2xs2gene_nomulti;
run;

data feat_on_off;
  set event.features_w_annotations_nomulti;
  if feature_type = 'fusion' then delete;
  keep feature_id flag_feature_on;
run;

data feat_to_keep;
   set event.flagged_feature_short;
   where flag_feature_short=0;
   keep feature_id;
run;

proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=feat_to_keep;
   by feature_id;
proc sort data=feat_on_off;
  by feature_id;
run;

data feat_on_off_w_xs;
   merge feat2xs2gene (in=in1) feat_on_off (in=in2) feat_to_keep (in=in3);
   by feature_id;
   if in1 and in2 and in3;
run;

* Calc number features on per xs;
proc sort data=feat_on_off_w_xs nodup;
   by transcript_id feature_id;
proc means data=feat_on_off_w_xs noprint;
   by transcript_id;
   var flag_feature_on;
   output out=num_feat_on_per_gene sum=num_features_dtct;
run;

data xs_w_no_feat_on;
   set num_feat_on_per_gene;
   where num_features_dtct = 0;
   keep transcript_id;
run;


/* Merge transcripts to remove */

proc sort data=xs_w_no_feat_on;
  by transcript_id;
proc sort data=xs_gene_no_exp;
  by transcript_id;
run;

data xscript2drop;
  merge xs_gene_no_exp (in=in1) xs_w_no_feat_on (in=in2);
  by transcript_id;
  if in1 then flag_xscript_gene_not_exp=1;
  else flag_xscript_gene_not_exp=0;
  if in2 then flag_xscript_no_feat_dtct=1;
  else flag_xscript_no_feat_dtct=0;
run;

proc freq data=xscript2drop;
  tables flag_xscript_gene_not_exp*flag_xscript_no_feat_dtct;
run;

/*

  flag_xscript_gene_not_exp
            flag_xscript_no_feat_dtct

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   2391 |   2391
           |   0.00 |  13.08 |  13.08
           |   0.00 | 100.00 |
           |   0.00 |  13.20 |
  ---------+--------+--------+
         1 |    170 |  15722 |  15892
           |   0.93 |  85.99 |  86.92
           |   1.07 |  98.93 |
           | 100.00 |  86.80 |
  ---------+--------+--------+
  Total         170    18113    18283
               0.93    99.07   100.00




*/

proc sort data=feat2xs2gene;
  by transcript_id;
proc sort data=xscript2drop;
  by transcript_id;
run;

data feat2xs2gene_keep;
  merge feat2xs2gene (in=in1) xscript2drop (in=in2);
  by transcript_id;
  if in1 and in2 then delete;
  else if in1 then output;
run;

/* Make permenant */

data event.feature2xs2gene_exp_only_nomulti;
  set feat2xs2gene_keep;
  drop flag_xscript_gene_not_exp flag_xscript_no_feat_dtct;
run;

