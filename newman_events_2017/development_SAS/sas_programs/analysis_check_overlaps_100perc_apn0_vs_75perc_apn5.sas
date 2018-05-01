ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* Compare overlap between 100% APN>0 list and 75% APN>5 list */

data xs_list1;
   set event.bin_xscripts_by_dtct_apn0;
   where perc_features_dtct=1;
   keep transcript_id;
run; *14734 transcripts;

data xs_list2;
   set event.bin_xscripts_by_dtct_apn5;
   where perc_features_dtct ge 0.75;
   keep transcript_id;
run; *13740 transcripts;

proc sort data=xs_list1;
   by transcript_id;
proc sort data=xs_list2;
   by transcript_id;
run;

data xs_overlap;
  merge xs_list1 (in=in1) xs_list2 (in=in2);
  by transcript_id;
  if in1 then flag_xscript_100perc_apn0=1; else flag_xscript_100perc_apn0=0;
  if in2 then flag_xscript_75perc_apn5=1; else flag_xscript_75perc_apn5=0;
run; *19164 transcripts total;

proc freq data=xs_overlap;
  tables flag_xscript_100perc_apn0*flag_xscript_75perc_apn5;
run;


/*
 flag_xscript_100perc_apn0
           flag_xscript_75perc_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |   4430 |   4430
          |   0.00 |  23.12 |  23.12
          |   0.00 | 100.00 |
          |   0.00 |  32.24 |
 ---------+--------+--------+
        1 |   5424 |   9310 |  14734
          |  28.30 |  48.58 |  76.88
          |  36.81 |  63.19 |
          | 100.00 |  67.76 |
 ---------+--------+--------+
 Total        5424    13740    19164
             28.30    71.70   100.00

68% transcripts (9310) in the 75% APN>5 list are in the 100% aPN>0

*/

/* Compare to PB: where do they sit? */

data pb2xs;
  set event.pacbio_to_refseq_nomulti_v2;
  keep transcript_id;
run;

proc sort data=xs_overlap;
  by transcript_id;
proc sort data=pb2xs nodup;
  by transcript_id;
run;

data xs_overlap_w_pb;
  merge xs_overlap (in=in1) pb2xs (in=in2);
  by transcript_id;
  if in2 then flag_in_pacbio=1; else flag_in_pacbio=0;
run;

proc freq data=xs_overlap_w_pb noprint;
  tables flag_xscript_75perc_apn5*flag_xscript_100perc_apn0*flag_in_pacbio / out=xs_count;
run;

proc print data=xs_count;
run;

/*
 flag_xscript_    flag_xscript_    flag_in_
  75perc_apn5      100perc_apn0     pacbio     COUNT

       .                .              1         874
       0                1              0        4820
       0                1              1         604
       1                0              0        3962
       1                0              1         468
       1                1              0        5376
       1                1              1        3934

3934 PB in both
604 in 100perc		77perc
468 in 75pec		


5006
*/


