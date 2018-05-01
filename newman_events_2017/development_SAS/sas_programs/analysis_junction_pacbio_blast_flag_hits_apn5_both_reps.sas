ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Flag hits with gaps or mismatch */

data flag_partial;
  set event.blast_junc_pb_hits_with_annot;
  if perc_identity >= 95 then flag_perc_ident_ge95=1; else flag_perc_ident_ge95=0;
  if mismatch = 0 then flag_no_mismatch=1; else flag_no_mismatch=0;
  if gapopen = 0 then flag_no_gap=1; else flag_no_gap=0;
run;

data dtct_events;
   set event.flag_feature_on_all_reps_ge5;
   keep feature_id flag_feature_on_ge5 flag_feature_on_ge10;
   rename feature_id=event_id;
run;

proc sort data=flag_partial;
  by event_id;
proc sort data=dtct_events;
  by event_id;
run;

data flag_partial2;
  merge flag_partial (in=in1) dtct_events (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=flag_partial2 noprint;
  tables flaG_perc_ident_ge95*flag_no_mismatch*flag_no_gap*flag_feature_on_ge5 / out=check;
proc print data=check;
run;

/*
flag_perc_    flag_no_    flag_no_    flag_feature_
ident_ge95    mismatch       gap          on_ge5        COUNT    PERCENT

     0            0           0             .              18      .
     0            0           0             0             166     0.1180
     0            0           0             1             171     0.1216
     0            0           1             .              18      .
     0            0           1             0             107     0.0761
     0            0           1             1             301     0.2140
     0            1           0             .               2      .
     0            1           0             0               7     0.0050
     0            1           0             1               9     0.0064
     1            0           0             .              50      .
     1            0           0             0             421     0.2994
     1            0           0             1             217     0.1543
     1            0           1             .             362      .
     1            0           1             0            2200     1.5644
     1            0           1             1            2166     1.5402
     1            1           0             .              37      .
     1            1           0             0             588     0.4181
     1            1           0             1             152     0.1081
     1            1           1             .           14720      .
     1            1           1             0           32900    23.3952
     1            1           1             1          101222    71.9791

101222 to keep
*/


data event_means;
   set event.mean_apn_events_npc;
   if mean_apn_npc > 0 then flag_mean_apn_gt0=1; else flag_mean_apn_gt0=0;
   if mean_apn_npc >= 5 then flag_mean_apn_ge5=1; else flag_mean_apn_ge5=0;
   if mean_apn_npc >= 10 then flag_mean_apn_ge10=1; else flag_mean_apn_ge10=0;
   rename feature_id=event_id;
run;

data event_annot;
   set event.features_w_annotations_nomulti;
   where feature_type="splicing" ;
   keep feature_id flag_feature_on flag_junction_annotated flag_intron_retention;
   rename feature_id=event_id;
run;

proc sort data=flag_partial2;
   by event_id;
proc sort data=event_means;
   by event_id;
proc sort data=event_annot;
   by event_id;
run;

data blast_hits_check_on;
  merge flag_partial2 (in=in1) event_annot (in=in2) event_means;
  by event_id;
  if in1 then flag_event_has_blast=1; else flag_event_has_blast=0;
  if in2 ;
run;


/* add fragment length */

data junc_length;
   set evspl.splicing_events_annot_refseq;
   feature_length=event_size;
   keep event_id feature_length gene_id;
run;

proc sort data=blast_hits_check_on;
  by event_id;
proc sort data=junc_length;
  by event_id;
run;

data blast_hits_w_len;
  merge blast_hits_check_on (in=in1) junc_length (in=in2);
  by event_id;
  if in2 ;
run;

data gene_exp;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

proc sort data=blast_hits_w_len;
  by gene_id;
proc sort data=gene_exp;
  by gene_id;
run;

data blast_hits_w_len2;
  merge blast_hits_w_len (in=in1) gene_exp (in=in2);
  by gene_id;
  if in1 and in2;
run;

data flag_length;
   set blast_hits_W_len2;
   if flag_event_has_blast=1 then do;
      if length >= 0.9 * feature_length then flag_hit_90perc_len=1; else flag_hit_90perc_len=0;
      if length >= 0.95 * feature_length then flag_hit_95perc_len=1; else flag_hit_95perc_len=0;
      if length >= 0.99 * feature_length then flag_hit_99perc_len=1; else flag_hit_99perc_len=0;
      if length >= 1 * feature_length then flag_hit_100perc_len=1; else flag_hit_100perc_len=0;
   end;
run;

*add these back in;
data event_dtct;
   set event.flag_feature_on_all_reps_ge5;
   keep feature_id flag_feature_on_ge5;
   rename feature_id=event_id;
run;

proc sort data=flag_length;
   by event_id;
proc sort data=event_dtct;
  by event_id;
run;

 data flag_length2;
   merge flag_length (in=in1) event_dtct (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=flag_length2;
   by event_id;
proc means data=flag_length2 noprint;
   by event_id;
   var flag_event_has_blast flag_pb_gene_has_refseq flag_junction_annotated flag_intron_retention
         flag_perc_ident_ge95 flag_no_mismatch flag_no_gap 
         flag_feature_on flag_hit_90perc_len flag_mean_apn_gt0 
         flag_mean_apn_ge5 flag_feature_on_ge5;
   output out=flag_length3 max=;
run;

/* count -- this is per hit */

proc freq data=flag_length2 noprint;
  tables flag_event_has_blast*flag_pb_gene_has_refseq*flag_junction_annotated*flag_intron_retention*
         flag_perc_ident_ge95*flag_no_mismatch*flag_no_gap*
         flag_feature_on*flag_hit_90perc_len*flag_mean_apn_gt0*
         flag_mean_apn_ge5*flag_feature_on_ge5 /
         out=pb_hit_count;
run;

proc export data=pb_hit_count outfile="!MCLAB/event_analysis/analysis_output/blast_junctions_to_pacbio_summary_of_hits.csv"
   dbms=csv replace;
run;

/* Extract only hits that are "good" */

data good_blast ;
  set flag_length3;
  if flag_pb_gene_has_refseq=1 and flag_perc_ident_ge95=1 and flag_no_mismatch=1 and flag_no_gap=1
  and flag_feature_on=1 and flag_hit_90perc_len=1 and flag_feature_on_ge5=1 ;
  keep event_id;
run;

data all_juncs;
   set flag_length3;
   if flag_feature_on_ge5=1;
   keep event_id flag_mean_apn_gt0 flag_mean_apn_ge5 flag_junction_annotated flag_intron_retention;
run;

proc sort data=good_blast nodup;
   by event_id;
proc sort data=all_juncs nodup;
  by event_id;
run;


data blast_summary;
  merge all_juncs (in=in1) good_blast (in=in2);
  by event_id;
  if in2 then flag_event_has_blast=1;
  else flag_event_has_blast=0;
run;

proc freq data=blast_summary noprint;
  tables  flag_mean_apn_gt0*flag_mean_apn_ge5*flag_junction_annotated*
          flag_intron_retention*flag_event_has_blast / out=summarize_flags;
run;

proc print data=summarize_flags;
run;

/*
  flag_      flag_
  mean_      mean_     flag_junction_    flag_intron_    flag_event_
 apn_gt0    apn_ge5       annotated        retention      has_blast     COUNT

    1          1              0                0              0             4
    1          1              0                0              1            16
    1          1              0                1              0           118
    1          1              0                1              1           214
    1          1              1                0              0           705
    1          1              1                0              1         27667

16/(16+4) unannotated with hits
214/(214+118) IR with hits
27667/(27667+705) annot with hits

*/

/* Make permanent */

data event.blast_junc_to_pb_summary;
  set blast_summary;
run;


data event.blast_junc_to_pb_all_hits;
  set flag_length;
run;


/* hits */

data hits;
  set flag_length2;
  if flag_perc_ident_ge95=1 and flag_no_mismatch=1
  and flag_no_gap=1 and	flag_hit_90perc_len=1;
  keep event_id;
run;

data annot;
  set evspl.splicing_events_annot_refseq;
  keep event_id gene_id flag_junction_annotated flag_intron_retention;
run;

data event_on;
   set event.flag_feature_on_all_reps_ge5;
   keep feature_id flag_feature_on_ge5;
   rename feature_id=event_id;
run;

data exp_Gene;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

proc sort data=exp_gene;
  by gene_id;
proc sort data=annot;
  by gene_id;
run;

data annot_exp;
  merge annot (in=in1) exp_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=event_on;
  by event_id;
proc sort data=annot_exp;
  by event_id;
proc sort data=hits nodup;
  by event_id;
run;

data hits_w_annot;
  merge annot_exp (in=in1) hits (in=in2) event_on (in=in3);
  by event_id;
  if in2 then flag_event_has_hit=1; else flag_event_has_hit=0;
  if in1;
run;

proc freq data=hits_w_annot noprint;
   tables flag_junction_annotated*flag_intron_retention*flag_feature_on_ge5*flag_event_has_hit / out=blast_check;
run;

proc print data=blast_check;
run;


/*
                                                     flag_
 flag_junction_    flag_intron_    flag_feature_     event_
    annotated        retention         on_ge5       has_hit      COUNT    PERCENT

        0                0               1             0           211     0.0112
        0                0               1             1           432     0.0229

        0                1               1             0           743     0.0394
        0                1               1             1           223     0.0118

        1                0               1             0          4627     0.2456
        1                0               1             1         27890     1.4804


        0                0               .             0           208      .
        0                0               .             1           148      .
        0                1               .             0           709      .
        0                1               .             1           167      .
        1                0               .             0          3258      .
        1                0               .             1          5734      .
        0                0               0             0       1563420    82.9837
        0                0               0             1           521     0.0277
        0                1               0             0        152495     8.0942
        0                1               0             1           650     0.0345
        1                0               0             0        121979     6.4744
        1                0               0             1         10818     0.5742

*/


data event.blast_junc_to_pb_w_annot;
  set hits_w_annot;
run;

