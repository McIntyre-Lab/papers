/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the simulation of 10000 transcripts, compare the list of selected transcripts
   against the transcript detection bins for APN>0 and APN>5
*/

data xs10000;
   set event.polyester_xs_list_60genes;
run;

data ea_apn0;
   set event.bin_xs_by_dtct_apn0_10gn;
   keep transcript_id bin_xscript_perc_dtct;
   rename bin_xscript_perc_dtct=bin_dtct_apn0;
run;

data ea_apn5;
   set event.bin_xs_by_dtct_apn5_10gn;
   keep transcript_id bin_xscript_perc_dtct;
   rename bin_xscript_perc_dtct=bin_dtct_apn5;
run;

proc sort data=xs10000;
  by transcript_id;
proc sort data=ea_apn0;
  by transcript_id;
proc sort data=ea_apn5;
  by transcript_id;
run;

data xs10000_v_ea;
   merge ea_apn0 (in=in1) ea_apn5 (in=in2) xs10000 (in=in3);
   by transcript_id;
   if in3 then flag_xscript_simulated=1;
   else flag_xscript_simulated=0;
   if in1 and in2 then output;
run;

proc freq data=xs10000_v_ea ;
  tables flag_xscript_simulated* bin_dtct_apn0  / out=sim_vs_ea_apn0;
  tables flag_xscript_simulated* bin_dtct_apn5  / out= sim_vs_ea_apn5;
run;

proc print data=sim_vs_ea_apn0;
run;

proc print data=sim_vs_ea_apn5;
run;

/*
flag_xscript_simulated     bin_dtct_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |1-24%   |100%    |50-74%  |75-99%  |  Total
---------+--------+--------+--------+--------+
       0 |      2 |      1 |      0 |      0 |      3
         |   0.43 |   0.21 |   0.00 |   0.00 |   0.64
         |  66.67 |  33.33 |   0.00 |   0.00 |
         | 100.00 |   0.23 |   0.00 |   0.00 |
---------+--------+--------+--------+--------+
       1 |      0 |    428 |      1 |     38 |    467
         |   0.00 |  91.06 |   0.21 |   8.09 |  99.36
         |   0.00 |  91.65 |   0.21 |   8.14 |
         |   0.00 |  99.77 | 100.00 | 100.00 |
---------+--------+--------+--------+--------+
Total           2      429        1       38      470
             0.43    91.28     0.21     8.09   100.00


flag_xscript_simulated     bin_dtct_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |0%      |1-24%   |100%    |50-74%  |75-99%  |  Total
---------+--------+--------+--------+--------+--------+
       0 |      1 |      2 |      0 |      0 |      0 |      3
         |   0.21 |   0.43 |   0.00 |   0.00 |   0.00 |   0.64
         |  33.33 |  66.67 |   0.00 |   0.00 |   0.00 |
         | 100.00 | 100.00 |   0.00 |   0.00 |   0.00 |
---------+--------+--------+--------+--------+--------+
       1 |      0 |      0 |    426 |      1 |     40 |    467
         |   0.00 |   0.00 |  90.64 |   0.21 |   8.51 |  99.36
         |   0.00 |   0.00 |  91.22 |   0.21 |   8.57 |
         |   0.00 |   0.00 | 100.00 | 100.00 | 100.00 |
---------+--------+--------+--------+--------+--------+
Total           1        2      426        1       40      470
             0.21     0.43    90.64     0.21     8.51   100.00


*/


