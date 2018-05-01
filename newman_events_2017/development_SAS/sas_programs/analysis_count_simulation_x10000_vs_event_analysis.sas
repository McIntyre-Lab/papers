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
   set event.polyester_xs_list_10k;
run;

data ea_apn0;
   set event.bin_xs_by_dtct_apn0_10k_v2;
   keep transcript_id bin_xscript_perc_dtct;
   rename bin_xscript_perc_dtct=bin_dtct_apn0;
run;

data ea_apn5;
   set event.bin_xs_by_dtct_apn5_10k_v2;
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
 Col Pct  |0%      |1-24%   |100%    |25-50%  |50-74%  |75-99%  |  Total
 ---------+--------+--------+--------+--------+--------+--------+
        0 |  46901 |   1082 |   5458 |   2043 |   6424 |  19910 |  81818
          |  51.08 |   1.18 |   5.94 |   2.23 |   7.00 |  21.68 |  89.11
          |  57.32 |   1.32 |   6.67 |   2.50 |   7.85 |  24.33 |
          | 100.00 | 100.00 |  37.05 |  99.61 |  99.15 |  96.77 |
 ---------+--------+--------+--------+--------+--------+--------+
        1 |      0 |      0 |   9272 |      8 |     55 |    665 |  10000
          |   0.00 |   0.00 |  10.10 |   0.01 |   0.06 |   0.72 |  10.89
          |   0.00 |   0.00 |  92.72 |   0.08 |   0.55 |   6.65 |
          |   0.00 |   0.00 |  62.95 |   0.39 |   0.85 |   3.23 |
 ---------+--------+--------+--------+--------+--------+--------+
 Total       46901     1082    14730     2051     6479    20575    91818
             51.08     1.18    16.04     2.23     7.06    22.41   100.00



  flag_xscript_simulated     bin_dtct_apn5

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |0%      |1-24%   |100%    |25-50%  |50-74%  |75-99%  |  Total
  ---------+--------+--------+--------+--------+--------+--------+
         0 |  47804 |    881 |   4380 |   2003 |   6929 |  19821 |  81818
           |  52.06 |   0.96 |   4.77 |   2.18 |   7.55 |  21.59 |  89.11
           |  58.43 |   1.08 |   5.35 |   2.45 |   8.47 |  24.23 |
           | 100.00 | 100.00 |  32.63 |  98.82 |  98.65 |  95.94 |
  ---------+--------+--------+--------+--------+--------+--------+
         1 |      0 |      0 |   9043 |     24 |     95 |    838 |  10000
           |   0.00 |   0.00 |   9.85 |   0.03 |   0.10 |   0.91 |  10.89
           |   0.00 |   0.00 |  90.43 |   0.24 |   0.95 |   8.38 |
           |   0.00 |   0.00 |  67.37 |   1.18 |   1.35 |   4.06 |
  ---------+--------+--------+--------+--------+--------+--------+
  Total       47804      881    13423     2027     7024    20659    91818
              52.06     0.96    14.62     2.21     7.65    22.50   100.00



*/


