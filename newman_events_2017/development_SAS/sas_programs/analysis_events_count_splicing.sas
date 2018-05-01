ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count the number of splicing that are expressed at APN>0 */

data splicing;
   set refseq.rfsq_flag_splicing_on;
   keep event_id flag_splicing_on;
run;

data annot;
   set evspl.splicing_events_annot_refseq;
   if transcript_id ne ' ' then flag_junction_annotated=1;
   keep event_id flag_junction_annotated flag_intron_retention
                 flag_exonskip flag_alt_donor flag_alt_acceptor;
   run;

proc sort data=splicing;
   by event_id;
proc sort data=annot;
   by event_id;
run;

data splicing_w_annot;
  merge splicing (in=in1) annot (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=splicing_w_annot;
   tables flag_splicing_on;
run;

/*

                                              Cumulative    Cumulative
 flag_splicing_on    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0     2818801       94.01       2818801        94.01
                1      179653        5.99       2998454       100.00
*/


/* Number of events detected by type */

proc freq data=splicing_w_annot;
  where flag_splicing_on=1;
  tables flag_junction_annotated*flag_intron_retention;
run;

/*
  flag_junction_annotated
            flag_intron_retention

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  12275 |  40443 |  52718
           |   6.83 |  22.51 |  29.34
           |  23.28 |  76.72 |
           |   8.82 | 100.00 |
  ---------+--------+--------+
         1 | 126935 |      0 | 126935
           |  70.66 |   0.00 |  70.66
           | 100.00 |   0.00 |
           |  91.18 |   0.00 |
  ---------+--------+--------+
  Total      139210    40443   179653
              77.49    22.51   100.00
*/

proc freq data=splicing_w_annot;
  where flag_splicing_on=1 and flag_intron_retention=0;
  tables flag_exonskip flag_alt_donor*flag_alt_acceptor;
run;

/*

                                             Cumulative    Cumulative
   flag_exonskip    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
               0      115182       82.74        115182        82.74
               1       24028       17.26        139210       100.00

  flag_alt_donor
            flag_alt_acceptor

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 127048 |   5254 | 132302
           |  91.26 |   3.77 |  95.04
           |  96.03 |   3.97 |
           |  96.20 |  73.54 |
  ---------+--------+--------+
         1 |   5018 |   1890 |   6908
           |   3.60 |   1.36 |   4.96
           |  72.64 |  27.36 |
           |   3.80 |  26.46 |
  ---------+--------+--------+
  Total      132066     7144   139210
              94.87     5.13   100.00

*/


/* Make permenant */

data event.splicing_on_apn_gt0;
   set splicing_w_annot;
run;

