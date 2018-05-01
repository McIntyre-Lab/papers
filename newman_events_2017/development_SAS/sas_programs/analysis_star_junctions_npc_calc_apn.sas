ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Convert STAR junction read counts to APN
   When running STAR, I had set the parameter '--sjdbOverhang' to 40
   this parameter is defined the length of the donor/acceptor sequence on each side of the junction
   This is the amount of overhang on either side So the total length of the junction sequence
   is 40*2=80, or about the same size as junctions in the catalog
 */

data star_juncs;
  set event.npc_star_junctions;
  drop flag_: ;
run;

data star_junc_apn;
   set star_juncs;
   NPC1_apn=NSC1/80;
   NPC2_apn=NSC2/80;
   mean_NPC_apn=(NPC1_apn+NPC2_apn)/2;
run;

data flag_apn;
   set star_junc_apn;
   if NPC1_apn > 0 then flag_NPC1_apn_gt0=1; else flag_NPC1_apn_gt0=0;
   if NPC2_apn > 0 then flag_NPC2_apn_gt0=1; else flag_NPC2_apn_gt0=0;
   if NPC1_apn ge 2 then flag_NPC1_apn_ge2=1; else flag_NPC1_apn_ge2=0;
   if NPC2_apn ge 2 then flag_NPC2_apn_ge2=1; else flag_NPC2_apn_ge2=0;
   if NPC1_apn ge 5 then flag_NPC1_apn_ge5=1; else flag_NPC1_apn_ge5=0;
   if NPC2_apn ge 5 then flag_NPC2_apn_ge5=1; else flag_NPC2_apn_ge5=0;
   if NPC1_apn ge 10 then flag_NPC1_apn_ge10=1; else flag_NPC1_apn_ge10=0;
   if NPC2_apn ge 10 then flag_NPC2_apn_ge10=1; else flag_NPC2_apn_ge10=0;
   if flag_NPC1_apn_gt0=1 and flag_NPC2_apn_gt0=1 then flaG_NPC_both_apn_gt0=1; else flaG_NPC_both_apn_gt0=0;
   if flag_NPC1_apn_ge2=1 and flag_NPC2_apn_ge2=1 then flaG_NPC_both_apn_ge2=1; else flaG_NPC_both_apn_ge2=0;
   if flag_NPC1_apn_ge5=1 and flag_NPC2_apn_ge5=1 then flaG_NPC_both_apn_ge5=1; else flaG_NPC_both_apn_ge5=0;
   if flag_NPC1_apn_ge10=1 and flag_NPC2_apn_ge10=1 then flaG_NPC_both_apn_ge10=1; else flaG_NPC_both_apn_ge10=0;
run;


data event.npc_star_junctions_est_apn;
   set flag_apn;
run;


proc freq data=flag_apn;
  tables flag_NPC1_apn_gt0*flag_NPC2_apn_gt0
         flag_NPC1_apn_ge2*flag_NPC2_apn_ge2
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
          0 |      0 |  17671 |  17671
            |   0.00 |  13.84 |  13.84
            |   0.00 | 100.00 |
            |   0.00 |  18.29 |
   ---------+--------+--------+
          1 |  31024 |  78967 | 109991
            |  24.30 |  61.86 |  86.16
            |  28.21 |  71.79 |
            | 100.00 |  81.71 |
   ---------+--------+--------+
   Total       31024    96638   127662
               24.30    75.70   100.00

flag_NPC1_apn_ge2
          flag_NPC2_apn_ge2

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 123915 |     10 | 123925
         |  97.06 |   0.01 |  97.07
         |  99.99 |   0.01 |
         |  98.70 |   0.47 |
---------+--------+--------+
       1 |   1635 |   2102 |   3737
         |   1.28 |   1.65 |   2.93
         |  43.75 |  56.25 |
         |   1.30 |  99.53 |
---------+--------+--------+
Total      125550     2112   127662
            98.35     1.65   100.00

 flag_NPC1_apn_ge5
           flag_NPC2_apn_ge5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 126500 |      0 | 126500
          |  99.09 |   0.00 |  99.09
          | 100.00 |   0.00 |
          |  99.53 |   0.00 |
 ---------+--------+--------+
        1 |    601 |    561 |   1162
          |   0.47 |   0.44 |   0.91
          |  51.72 |  48.28 |
          |   0.47 | 100.00 |
 ---------+--------+--------+
 Total      127101      561   127662
             99.56     0.44   100.00


flag_NPC1_apn_ge10
          flag_NPC2_apn_ge10

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 127226 |      0 | 127226
         |  99.66 |   0.00 |  99.66
         | 100.00 |   0.00 |
         |  99.80 |   0.00 |
---------+--------+--------+
       1 |    250 |    186 |    436
         |   0.20 |   0.15 |   0.34
         |  57.34 |  42.66 |
         |   0.20 | 100.00 |
---------+--------+--------+
Total      127476      186   127662
            99.85     0.15   100.00

*/
