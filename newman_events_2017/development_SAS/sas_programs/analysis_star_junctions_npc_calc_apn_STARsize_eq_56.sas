ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Convert STAR junction read counts to APN
   When running STAR, I had set the parameter '--sjdbOverhang' to 55
   this parameter is defined the length of the donor/acceptor sequence on each side of the junction
   This is the amount of overhang on either side, but since there is no indication
   about the actual junction size (or indeed how junctions are mapping),
   I will set the junction size to the read size (56bp), which is the minimum junction size
 */

data star_juncs;
  set event.npc_star_junctions;
  drop flag_: ;
run;

data star_junc_apn;
   set star_juncs;
   NPC1_apn=NSC1/56;
   NPC2_apn=NSC2/56;
   mean_NPC_apn=(NPC1_apn+NPC2_apn)/2;
run;

data flag_apn;
   set star_junc_apn;
   if NPC1_apn > 0 then flag_NPC1_apn_gt0=1; else flag_NPC1_apn_gt0=0;
   if NPC2_apn > 0 then flag_NPC2_apn_gt0=1; else flag_NPC2_apn_gt0=0;
   if NPC1_apn ge 5 then flag_NPC1_apn_ge5=1; else flag_NPC1_apn_ge5=0;
   if NPC2_apn ge 5 then flag_NPC2_apn_ge5=1; else flag_NPC2_apn_ge5=0;
   if NPC1_apn ge 10 then flag_NPC1_apn_ge10=1; else flag_NPC1_apn_ge10=0;
   if NPC2_apn ge 10 then flag_NPC2_apn_ge10=1; else flag_NPC2_apn_ge10=0;
   if flag_NPC1_apn_gt0=1 and flag_NPC2_apn_gt0=1 then flaG_NPC_both_apn_gt0=1; else flaG_NPC_both_apn_gt0=0;
   if flag_NPC1_apn_ge5=1 and flag_NPC2_apn_ge5=1 then flaG_NPC_both_apn_ge5=1; else flaG_NPC_both_apn_ge5=0;
   if flag_NPC1_apn_ge10=1 and flag_NPC2_apn_ge10=1 then flaG_NPC_both_apn_ge10=1; else flaG_NPC_both_apn_ge10=0;
run;


data event.npc_star_junctions_est_apn_56bp;
   set flag_apn;
run;


proc freq data=flag_apn;
  tables flag_NPC1_apn_gt0*flag_NPC2_apn_gt0
         flag_NPC1_apn_ge5*flag_NPC2_apn_ge5
         flag_NPC1_apn_ge10*flag_NPC2_apn_ge10;
run;

/*
   flag_NPC1_apn_gt0
             flag_NPC2_apn_gt0

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  19111 |  23948 |  43059
            |   9.80 |  12.28 |  22.08
            |  44.38 |  55.62 |
            |  30.99 |  17.96 |
   ---------+--------+--------+
          1 |  42564 | 109417 | 151981
            |  21.82 |  56.10 |  77.92
            |  28.01 |  71.99 |
            |  69.01 |  82.04 |
   ---------+--------+--------+
   Total       61675   133365   195040
               31.62    68.38   100.00

 flag_NPC1_apn_ge5
           flag_NPC2_apn_ge5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 187774 |     10 | 187784
          |  96.27 |   0.01 |  96.28
          |  99.99 |   0.01 |
          |  98.38 |   0.24 |
 ---------+--------+--------+
        1 |   3091 |   4165 |   7256
          |   1.58 |   2.14 |   3.72
          |  42.60 |  57.40 |
          |   1.62 |  99.76 |
 ---------+--------+--------+
 Total      190865     4175   195040
             97.86     2.14   100.00

flag_NPC1_apn_ge10
          flag_NPC2_apn_ge10

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 191755 |      4 | 191759
         |  98.32 |   0.00 |  98.32
         | 100.00 |   0.00 |
         |  99.19 |   0.23 |
---------+--------+--------+
       1 |   1567 |   1714 |   3281
         |   0.80 |   0.88 |   1.68
         |  47.76 |  52.24 |
         |   0.81 |  99.77 |
---------+--------+--------+
Total      193322     1718   195040
            99.12     0.88   100.00

*/
