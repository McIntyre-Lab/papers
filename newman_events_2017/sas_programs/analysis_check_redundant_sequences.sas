
libname event '!MCLAB/event_analysis/sas_data';

*check: how many redundant sequences have coverage?;

data events_on;
  set event.splicing_on_apn_gt0;
  keep event_id flag_splicing_on;
run;


data seq_flag;
  set splice.flag_event_redundant_seq;
run;

proc sort data=events_on;
   by event_id;
proc sort data=seq_flag;
   by event_id;
run;

data events_on_seq_flag;
  merge events_on (in=in1) seq_flag (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=events_on_seq_flag;
   tables flag_splicing_on*flag_redundant_sequence;
run;


/*
    flag_splicing_on
              flag_redundant_sequence

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |2742852 |  75949 |2818801
             |  91.48 |   2.53 |  94.01
             |  97.31 |   2.69 |
             |  93.85 | 100.00 |
    ---------+--------+--------+
           1 | 179653 |      0 | 179653
             |   5.99 |   0.00 |   5.99
             | 100.00 |   0.00 |
             |   6.15 |   0.00 |
    ---------+--------+--------+
    Total     2922505    75949  2998454
                97.47     2.53   100.00


Okay, so all redundant sequences have no coverage, as thought.
*/
