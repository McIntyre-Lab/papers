ods listing; ods html close;

libname eventloc '/mnt/store/event_sandbox/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* Combine all splicing on flags for case-only data */

data all_flags;
  set eventloc.flag_splicing_ge5_: ;
run;

/* Flag if on in all cell types */

data flag_on_all;
  set all_flags;
  if flag_cd19_on=1 and flag_cd4_on=1 and flag_cd8_on=1 then flag_splicing_all_ge5=1;
  else flag_splicing_all_ge5=0;
run;

proc freq data=flag_on_all;
  tables flag_cd19_on flag_cd4_on flag_cd8_on flag_splicing_all_ge5;
run;

/*

                                            Cumulative    Cumulative
   flag_CD19_on    Frequency     Percent     Frequency      Percent
   -----------------------------------------------------------------
              0      377344       78.62        377344        78.62
              1      102592       21.38        479936       100.00

                       Frequency Missing = 236395


                                            Cumulative    Cumulative
    flag_CD4_on    Frequency     Percent     Frequency      Percent
    ----------------------------------------------------------------
              0      416305       80.05        416305        80.05
              1      103738       19.95        520043       100.00

                       Frequency Missing = 196288

      flag_CD8_on    Frequency     Percent     Frequency      Percent
      ----------------------------------------------------------------
                0      426726       80.16        426726        80.16
                1      105624       19.84        532350       100.00

                         Frequency Missing = 183981


     flag_splicing_                             Cumulative    Cumulative
            all_ge5    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0      622635       86.92        622635        86.92
                  1       93696       13.08        716331       100.00

*/

proc freq data=flag_on_all noprint;
  tables flag_cd19_on*flag_cd4_on*flag_cd8_on / out=junc_count;
proc print data=junc_count;
run;

/*

 flag_      flag_     flag_
CD19_on    CD4_on    CD8_on     COUNT

   .          .         .       39982
   .          .         0       62936
   .          .         1        2205
   .          0         .       53811
   .          0         0       68603
   .          0         1          35
   .          1         .        1209
   .          1         0           3
   .          1         1        7611
   0          .         .       50604
   0          .         0       33684
   0          .         1         217
   0          0         .       32081
   0          0         0      260327
   0          0         1          35
   0          1         .         188
   0          1         0           9
   0          1         1         199
   1          .         .        4850
   1          .         0         188
   1          .         1        1622
   1          0         .         433
   1          0         0         976
   1          0         1           4
   1          1         .         823
   1          1         1       93696

*/

/* Make permenant */

data event.t1d_flag_splicing_on_apn5;
  set flag_on_all;
  keep event_id flag_cd19_on flag_cd4_on flag_cd8_on flag_splicing_all_ge5;
run;

