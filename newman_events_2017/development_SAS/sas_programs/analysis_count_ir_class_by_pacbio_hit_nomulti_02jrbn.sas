ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count reclassified IR events by type and by whether they have a PB hit or not */

proc freq data=event.ir_reclassification_v3_apn5 noprint;
  tables flag_low_expressed*flag_possible_unprocessed*flag_possible_novel_donor*flag_possible_ir / out=check;
run;

proc print data=check;
run;



* IR with hit;
data ir_w_hit;
  set event.blast_dtct_jnc2pb_nomult_summary;
  where flag_correct_gene=1 and flag_feature_on=1 and flag_intron_retention=1;
  keep event_id;
run;

* IR reclassification;
data ir_class;
  set event.ir_reclassification_v3_apn5;
  where flag_splicing_on=1 and flag_low_expressed=0;
  keep event_id flag_low_expressed flag_possible_unprocessed flag_possible_novel_donor flag_possible_ir ;
run;


proc sort data=ir_w_hit nodup;
   by event_id;
proc sort data=ir_class;
   by event_id;
run;

data ir_class_pb_hit;
  merge ir_class (in=in1) ir_w_hit (in=in2);
  if in2 then flag_pb_hit=1; else flag_pb_hit=0;
  if in1 then output;
run;

proc freq data=ir_class_pb_hit noprint;
   tables flag_low_expressed*flag_possible_unprocessed*
          flag_possible_novel_donor*flag_possible_ir*flag_pb_hit / out=ir_count;
run;

proc print data=ir_count;
run;

/*
                                                      flag_

                                                   flag_
flag_low_    flag_possible_    flag_possible_    possible_    flag_pb_
expressed      unprocessed       novel_donor         ir          hit      COUNT    PERCENT

    0               0                 0              1            0         194     0.7652
    0               0                 0              1            1           6     0.0237
    0               0                 1              1            0          42     0.1657
    0               0                 1              1            1           2     0.0079
    0               1                 0              0            0       24208    95.4875
    0               1                 0              0            1         586     2.3115
    0               1                 1              0            0         282     1.1123
    0               1                 1              0            1          32     0.1262
    1               .                 .              .            0        1491      .
    1               .                 .              .            1         178      .

TABLE

CLASS		NO HIT	HIT
Poss. IR	194	6
Ambig IR	42	2
Unprocessed	24208	586
Novel donor	282	32
Low Express	1491	178

*/


