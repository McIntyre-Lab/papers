Total number of fusions: 275177
Number of fusions shorter than read length: 17591 (6%)
Number of fusions shorter than half of read length: 3592 (1%)


Total number of exons: 355269
Number of exons shorter than read length: 23027 (6%)
Number of exons shorter than half of read length: 4825 (1%)


Total number of exon fragments: 360236
Number of exon fragments shorter than read length: 64649 (17%)
Number of exon fragments shorter than half of read length: 38172 (10%)
Number of exon fragments 1bp in length: 5595 (1%)



/* Detection of fusions by commonality type */

data fus_on_flags;
  set isoold.rfsq_flag_fusions_on;
  keep fusion_id flag_fusion_on;
run;


data fus_com fus_con fus_uniq;
  set isoold.ref_flag_feature_commonality;
  where feature_type="fusion";
  if feature_group="unique" then output fus_uniq;
  else if feature_group="common" then output fus_com;
  else if feature_group="constitutive" then output fus_con;
  keep feature_id;
  rename feature_id=fusion_id;
run;

proc sort data=fus_uniq nodup;
  by fusion_id;
proc sort data=fus_com  nodup;
  by fusion_id;
proc sort data=fus_con  nodup;
  by fusion_id;
run;

data fus_com2;
  merge fus_uniq (in=in1) fus_com (in=in2) fus_con (in=in3);
  length feature_type $20.;
  by fusion_id;
  if in1 then feature_type="unique";
  else if in2 then feature_type="common";
  else if in3 then feature_type="constitutive";
run;


proc sort data=fus_on_flags;
   by fusion_id;
proc sort data=fus_com2 nodup;
   by fusion_id;
run;

data fus_on_com;
  merge fus_on_flags (in=in1) fus_com2 (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc freq data=fus_on_com noprint;
  tables feature_type*flag_fusion_on / out=fus_on_sum;
run;

proc print data=fus_on_sum;
run;




/* Detection of junctions by commonality type */

data splice_on_flags;
  set isoold.rfsq_flag_splicing_on;
  keep event_id flag_splicing_on;
run;


data splice_com splice_con splice_uniq;
  set isoold.ref_flag_feature_commonality;
  where feature_type="junction";
  if feature_group="unique" then output splice_uniq;
  else if feature_group="common" then output splice_com;
  else if feature_group="constitutive" then output splice_con;
  keep feature_id;
  rename feature_id=event_id;
run;

proc sort data=splice_uniq nodup;
  by event_id;
proc sort data=splice_com  nodup;
  by event_id;
proc sort data=splice_con  nodup;
  by event_id;
run;

data splice_com2;
  merge splice_uniq (in=in1) splice_com (in=in2) splice_con (in=in3);
  length feature_type $20.;
  by event_id;
  if in1 then feature_type="unique";
  else if in2 then feature_type="common";
  else if in3 then feature_type="constitutive";
run;

proc sort data=splice_on_flags;
   by event_id;
proc sort data=splice_com2;
   by event_id;
run;

data splice_on_com;
  merge splice_on_flags (in=in1) splice_com2 (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=splice_on_com noprint;
  tables feature_type*flag_splicing_on / out=splice_on_sum;
run;

proc print data=splice_on_sum;
run;

/* Detection of novel events */

data novel;
  set splicing.splicing_events_annot_refseq;
  if num_transcripts > 0 then flag_junction_annotated=1;
  keep event_id gene_id flag_junction_annotated flag_intron_retention;
run;

proc sort data=novel;
   by event_id;
proc sort data=splice_on_flags;
  by event_id;
run;

data novel_on;
  merge novel (in=in1) splice_on_flags (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=novel_on noprint;
  tables flag_junction_annotated*flag_intron_retention*flag_splicing_on / out=novel_on_sum;
run;

proc print data=novel_on_sum;
run;


data on_fus;
  set isoold.rfsq_flag_fusions_on;
  if flag_fusion_on=1;
  keep fusion_id;
  rename fusion_id=feature_id;
run;


data on_splice;
  set isoold.rfsq_flag_splicing_on;
  if flag_splicing_on=1;
  keep event_id;
  rename event_id=feature_id;
run;


data on_feat;
  set on_splice on_fus;
run;

data feat2xs;
  set isoold.ref_fusions_junctions_stack;
run;

proc sort data=feat2xs;
  by feature_id;
proc sort data=on_feat;
  by feature_id;
run;

data xs_on;
  merge feat2xs (in=in1)  on_feat (in=in2);
  by feature_id;
  if in1 and in2;
  keep transcript_id gene_id;
run;


data xs_on2;
  set xs_on;
  keep transcript_id;
run;

proc sort data=xs_on2 nodup;
  by transcript_id;
run;


data gene_on2;
  set xs_on;
  keep gene_id;
run;

proc sort data=gene_on2 nodup;
  by gene_id;
run;


data xs_uniq;
  set isoold.ref_flag_xs_w_uniq_feature;
  where flag_isoform_has_uniq_feat=1;
  keep transcript_id;
run;


data gene_uniq;
  set isoold.ref_flag_xs_w_uniq_feature;
  where flag_isoform_has_uniq_feat=1;
  keep gene_id;
run;

proc sort data=xs_uniq nodup;
  by transcript_id;
proc sort data=xs_on2 nodup;
  by transcript_id;
proc sort data=gene_uniq nodup;
  by gene_id;
proc sort data=gene_on2 nodup;
  by gene_id;
run;

data xs_on_uniq;
  merge xs_uniq (in=in1) xs_on2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene_on_uniq;
  merge gene_uniq (in=in1) gene_on2 (in=in2);
  by gene_id;
  if in1 and in2;
run;

data uniq_feat;
  set isoold.ref_flag_feature_commonality;
  where feature_group="unique";
  keep feature_id;
run;

proc sort data=uniq_feat;
  by feature_id;
proc sort data=on_feat;
  by feature_id;
run;

data on_uniq;
  merge on_feat (in=in1) uniq_feat (in=in2);
  by feature_id;
  if in1 and in2;
run;

proc sort data=feat2xs;
  by feature_id;
proc sort data=on_uniq;
  by feature_id;
run;

data xs_w_on_uniq;
  merge feat2xs (in=in1) on_uniq (in=in2);
  by feature_id;
  if in1 and in2;
  keep transcript_id gene_id;
run;

data xs_on_uniq;
  set xs_w_on_uniq;
  keep transcript_id;
run;

data gene_on_uniq;
  set xs_w_on_uniq;
  keep gene_id;
run;


proc sort data=xs_on_uniq nodup;
  by transcript_id;
proc sort data=gene_on_uniq nodup;
  by gene_id;
run;



/*
                        flag_
                       fusion_
Obs    feature_type       on      COUNT    PERCENT

 1     common             0       25731     9.3507
 2     common             1       34632    12.5854
 3     constitutive       0       49895    18.1320
 4     constitutive       1       86346    31.3785
 5     unique             0       39130    14.2200
 6     unique             1       39442    14.3334

Total on: 160420
Total off: 114756
Total: 275176



                          flag_
        feature_        splicing_
 Obs    group               on       COUNT    PERCENT

  1     common              0        45131    16.2225
  2     common              1        33493    12.0392
  3     constitutive        0        47489    17.0701
  4     constitutive        1        57554    20.6880
  5     unique              0        62346    22.4105
  6     unique              1        32187    11.5697


Total on: 123234
Total off: 154966
Total: 278200


                                            flag_
        flag_junction_    flag_intron_    splicing_
 Obs       annotated        retention         on        COUNT     PERCENT

  1            0                0             0        2448141    81.6468
  2            0                0             1          12275     0.4094
  3            0                1             0         207540     6.9216
  4            0                1             1          40443     1.3488


"on" fusions and splicing events represent:
98987 transcripts,
26571 genes

# xs with uniq			55684
# genes with uniq		29348

# xs with uniq on		29649
# genes with uniq on		17456

*/

How many transcripts do these represent?



data novel;
  set splicing.splicing_events_annot_refseq;
  if num_transcripts > 0 then flag_junction_annotated=1;
  keep event_id gene_id flag_junction_annotated flag_intron_retention;
run;

proc sort data=novel;
   by event_id;
proc sort data=splice_on_flags;
  by event_id;
run;

data novel_on;
  merge novel (in=in1) splice_on_flags (in=in2);
  by event_id;
  if in1 and in2;
run;


data genes;
  set novel_on;
  where flag_splicing_on=1 and flag_junction_annotated=0;
  keep gene_id ;
run;

proc sort data=genes nodup;
  by gene_id;
run;

