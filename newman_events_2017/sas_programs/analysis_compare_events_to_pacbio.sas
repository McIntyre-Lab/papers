ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Comparing the detected splicing events (from the lists of filtered transcripts and genes) to that of junctions derived from PacBio transcripts

Also want to count the number of PB transcripts that correspond to the transcripts in each bin*/

/* Check transcripts in common */
* Import list of PacBio transcript;

    data WORK.PB_XSCRIPT    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/useful_mouse_data/mm10/gff/pacbio/pacbio_transcript_list.txt'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat transcript_id $18. ;
       format transcript_id $18. ;
    input
                transcript_id $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

data xscripts_by_bin;
  set event.xscripts_w_unique_by_bin;
  keep transcript_id bin_xscript_perc_uniq_dtct;
run;

proc sort data=pb_xscript nodup;
   by transcript_id;
proc sort data=xscripts_by_bin;
   by transcript_id;
run;

data xscript_bin_w_pb;
  merge xscripts_by_bin (in=in1) pb_xscript (in=in2);
  by transcript_id;
  if in1 and in2 then do;
      flag_pacbio_match=1;
      output; end;
  else if in1 then do;
    flag_pacbio_match=0;
      output; end;
run;

proc freq data=xscript_bin_w_pb;
   tables flag_pacbio_match;
run;

/*
                               The FREQ Procedure

         flag_pacbio_                             Cumulative    Cumulative
                match    Frequency     Percent     Frequency      Percent
     ---------------------------------------------------------------------
                    0       93577       91.99         93577        91.99
                    1        8143        8.01        101720       100.00

*/


proc freq data=xscript_bin_w_pb;
   tables flag_pacbio_match*bin_xscript_perc_uniq_dtct;
run;

/*
            Table of flag_pacbio_match by bin_xscript_perc_uniq_dtct

flag_pacbio_match     bin_xscript_perc_uniq_dtct

Frequency|
Percent  |
Row Pct  |
Col Pct  |0%      |0-25%   |100%    |25-50%  |50-75%  |75-100% |no uniqu|  Total
         |        |        |        |        |        |        |e       |
---------+--------+--------+--------+--------+--------+--------+--------+
       0 |  22905 |   1655 |  12438 |   2742 |   7769 |    798 |  45270 |  93577
         |  22.52 |   1.63 |  12.23 |   2.70 |   7.64 |   0.78 |  44.50 |  91.99
         |  24.48 |   1.77 |  13.29 |   2.93 |   8.30 |   0.85 |  48.38 |
         |  98.18 |  99.22 |  72.95 |  98.88 |  95.80 |  81.02 |  94.70 |
---------+--------+--------+--------+--------+--------+--------+--------+
       1 |    424 |     13 |   4612 |     31 |    341 |    187 |   2535 |   8143
         |   0.42 |   0.01 |   4.53 |   0.03 |   0.34 |   0.18 |   2.49 |   8.01
         |   5.21 |   0.16 |  56.64 |   0.38 |   4.19 |   2.30 |  31.13 |
         |   1.82 |   0.78 |  27.05 |   1.12 |   4.20 |  18.98 |   5.30 |
---------+--------+--------+--------+--------+--------+--------+--------+
Total       23329     1668    17050     2773     8110      985    47805   101720
            22.93     1.64    16.76     2.73     7.97     0.97    47.00   100.00

*/



/* Now compare junctions detected */

* Import PB junctions;



5. Compare PB junctions to catalog:
	- all events
	- detected events
	- probable events

