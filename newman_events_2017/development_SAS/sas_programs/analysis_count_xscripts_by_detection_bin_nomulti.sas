ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Need to get counts and percentages for a table:

XS w/ uniq, XS w/o uniq
 by
% dtct bin (1-25%, 25-50%, 50%-75%, 75%-99%, 100%)

Then, for each cell of the table, what proportion have a PB match ? */

data xs_by_dtct;
   set event.xscripts_w_unique_by_bin;
   length bin_xscript_perc_dtct $10.;
   if perc_features_dtct = 1 then bin_xscript_perc_dtct="100%";
   else if perc_features_dtct >= 0.75 then bin_xscript_perc_dtct="75-99%";
   else if perc_features_dtct >= 0.5 then bin_xscript_perc_dtct="50-74%";
   else if perc_features_dtct >= 0.25 then bin_xscript_perc_dtct="25-49%";
   else if perc_features_dtct > 0 then bin_xscript_perc_dtct="1-24%";
   else bin_xscript_perc_dtct="0%";
   keep transcript_id flag_xscript_has_unique bin_xscript_perc_dtct;
run;

proc freq data=xs_by_dtct ;
  tables flag_xscript_has_unique*bin_xscript_perc_dtct;
run;

/*
 flag_xscript_has_unique     bin_xscript_perc_dtct

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |1-24%   |100%    |25-49%  |50-74%  |75-99%  |  Total
 ---------+--------+--------+--------+--------+--------+
        0 |   9342 |   5004 |   6628 |   5837 |   8089 |  34900
          |  12.70 |   6.80 |   9.01 |   7.94 |  11.00 |  47.46
          |  26.77 |  14.34 |  18.99 |  16.72 |  23.18 |
          |  58.52 |  33.96 |  56.70 |  51.83 |  40.67 |
 ---------+--------+--------+--------+--------+--------+
        1 |   6621 |   9730 |   5061 |   5424 |  11799 |  38635
          |   9.00 |  13.23 |   6.88 |   7.38 |  16.05 |  52.54
          |  17.14 |  25.18 |  13.10 |  14.04 |  30.54 |
          |  41.48 |  66.04 |  43.30 |  48.17 |  59.33 |
 ---------+--------+--------+--------+--------+--------+
 Total       15963    14734    11689    11261    19888    73535
             21.71    20.04    15.90    15.31    27.05   100.00

*/

/* Now count number of transcripts with PB hits */

data xs_w_pb;
   set event.pacbio2refseq_id_nomulti;
   keep transcript_id;
run;

proc sort data=xs_w_pb nodup;
  by transcript_id;
proc sort data=xs_by_dtct;
  by transcript_id;
run;

data xs_by_dtct_w_pb;
  merge xs_by_dtct (in=in1) xs_w_pb (in=in2);
  by transcript_id;
  if in2 then flag_has_pacbio=1; else flag_has_pacbio=0;
  if in1;
run;

proc freq data=xs_by_dtct_w_pb;
  where flag_has_pacbio=1;
  tables flag_xscript_has_unique*bin_xscript_perc_dtct;
run;

/*
 flag_xscript_has_unique     bin_xscript_perc_dtct

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |1-24%   |100%    |25-49%  |50-74%  |75-99%  |  Total
 ---------+--------+--------+--------+--------+--------+
        0 |      4 |    953 |     10 |     21 |    350 |   1338
          |   0.09 |  21.31 |   0.22 |   0.47 |   7.82 |  29.91
          |   0.30 |  71.23 |   0.75 |   1.57 |  26.16 |
          |  21.05 |  28.87 |  41.67 |  17.21 |  34.76 |
 ---------+--------+--------+--------+--------+--------+
        1 |     15 |   2348 |     14 |    101 |    657 |   3135
          |   0.34 |  52.49 |   0.31 |   2.26 |  14.69 |  70.09
          |   0.48 |  74.90 |   0.45 |   3.22 |  20.96 |
          |  78.95 |  71.13 |  58.33 |  82.79 |  65.24 |
 ---------+--------+--------+--------+--------+--------+
 Total          19     3301       24      122     1007     4473
              0.42    73.80     0.54     2.73    22.51   100.00
*/

