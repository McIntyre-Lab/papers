ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Last check: I want to look at the cross-over between genes with unannotated hits and gets with novel PB transcripts */


data unannot_feat;
   set event.features_w_annotations_nomulti;
   where flag_feature_on=1 and flag_junction_annotated=0 and feature_type="splicing";
   keep feature_id flag_intron_retention;
run;

data feat2gene;
  set evspl.splicing_Events_annot_refseq;
  where num_transcripts=0;
  keep event_id gene_id;
  rename event_id=feature_id;
run;

data pb2gene;
   set event.pacbio2refseq_gene_nomulti;
   keep pacbio_gene_id gene_id;
run;

data pb_novel;
   set event.pacbio_transcripts_nomulti;
   where pb_status="Novel";
   keep pacbio_gene_id;
run;

proc sort data=unannot_feat;
  by feature_id;
proc sort data=feat2gene;
  by feature_id;
run;

data unannot2gene;
  merge feat2gene (in=in1) unannot_feat (in=in2);
  by feature_id;
  if in1 and in2;
run;

proc sort data=pb2gene nodup;
   by pacbio_gene_id gene_id;
proc sort data=pb_novel nodup;
   by pacbio_gene_id;
run;

data pb2gene_novel;
  merge pb2gene (in=in1) pb_novel (in=in2);
  by pacbio_gene_id;
  if in1 and in2;
  keep gene_id;
run;

data gene_w_uannot;
   set unannot2gene;
   keep gene_id;
run;


proc sort data=pb2gene_novel nodup;
   by gene_id;
proc sort data=gene_w_uannot nodup;
  by gene_id;
run;

data unannot_flag_pb_novel;
  merge gene_w_uannot (in=in1) pb2gene_novel (in=in2);
  by gene_id;
  if in1 then flag_has_unannot=1; else flag_has_unannot=0;
  if in2 then flag_has_pb_novel=1; else flag_has_pb_novel=0;
run;

proc freq data=unannot_flag_pb_novel;
  tables flag_has_unannot*flag_has_pb_novel;
run;

/*
 flag_has_unannot
           flag_has_pb_novel

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |     83 |     83
          |   0.00 |   0.83 |   0.83
          |   0.00 | 100.00 |
          |   0.00 |   4.64 |
 ---------+--------+--------+
        1 |   8204 |   1704 |   9908
          |  82.11 |  17.06 |  99.17
          |  82.80 |  17.20 |
          | 100.00 |  95.36 |
 ---------+--------+--------+
 Total        8204     1787     9991
             82.11    17.89   100.00

*/

data exp_genes;
  set event.flag_gene_expressed;
  where flaG_gene_expressed=1;
  keep gene_id;
run;

proc sort data=unannot_flag_pb_novel;
   by gene_id;
proc sort data=exp_genes nodup;
  by gene_id;
run;

data unannot_flag_pb_novel_exp;
  merge unannot_flag_pb_novel (in=in1) exp_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc freq data=unannot_flag_pb_novel_exp;
  tables flag_has_unannot*flag_has_pb_novel;
run;

/*
Table of flag_has_unannot by flag_has_pb_novel

     flag_has_unannot
               flag_has_pb_novel

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |      0 |     82 |     82
              |   0.00 |   0.86 |   0.86
              |   0.00 | 100.00 |
              |   0.00 |   4.59 |
     ---------+--------+--------+
            1 |   7737 |   1704 |   9441
              |  81.25 |  17.89 |  99.14
              |  81.95 |  18.05 |
              | 100.00 |  95.41 |
     ---------+--------+--------+
     Total        7737     1786     9523
                 81.25    18.75   100.00
*/


