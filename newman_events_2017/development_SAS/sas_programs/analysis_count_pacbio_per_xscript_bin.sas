ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* For the various transcript bins, how many have a matching Pacbio transcript? */

* transcript bins;

data xscript_bins;
  set event.xscripts_w_unique_by_bin;
  keep transcript_id bin_xscript_perc_uniq_dtct;
run;

* Get list of all transcripts;

data refseq_xscripts;
   set event.feature2xs2gene;
   keep transcript_id;
run;

* pacbio2refseq;
data pb2refseq;
  set event.pacbio2refseq_id;
run;


/* sort and merge */
proc sort data=xscript_bins;
  by transcript_id;
proc sort data=refseq_xscripts nodup;
  by transcript_id;
proc sort data=pb2refseq nodup;
  by transcript_id;
run;

data xscript_bins_all;
  merge xscript_bins (in=in1) refseq_xscripts (in=in2);
  by transcript_id;
  if not in1 then bin_xscript_perc_uniq_dtct="not_exp";
run;

data xscript_bin2pacbio;
  merge xscript_bins_all (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 and in2 then flag_xscript_has_pacbio=1;
  else if in1 then flag_xscript_has_pacbio=0;
  else do;
     flag_xscript_has_pacbio=.;
     bin_xscript_perc_uniq_dtct="no_refseq";
  end;
run;


/* Check: how many matches? */

proc freq data=xscript_bin2pacbio;
   tables flag_xscript_has_pacbio;
run;

/*
  flag_xscript_
     has_pacbio    Frequency     Percent
-----------------------------------------
              0      119785       92.54
              1        9653        7.46
*/

/* There appear to be some Refseq transcripts that match more than one PB transcript
   so I am only going to consider if the Refseq transcript has any PB, rather than a 1-to-1 match */

data xscript_check;
   set xscript_bin2pacbio;
   keep transcript_id  bin_xscript_perc_uniq_dtct flag_xscript_has_pacbio;
run;

proc sort data=xscript_check nodup;
  by transcript_id;
run;

proc freq data=xscript_check;
   tables flag_xscript_has_pacbio*bin_xscript_perc_uniq_dtct;
run;


/*
  flag_xscript_has_pacbio     bin_xscript_perc_uniq_dtct

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |0%      |0-25%   |100%    |25-50%  |50-75%  |75-100% |no uniqu|not_exp |  Total
           |        |        |        |        |        |        |e       |        |
  ---------+--------+--------+--------+--------+--------+--------+--------+--------+
         0 |  22844 |   1655 |  12215 |   2719 |   7687 |    732 |  45105 |  26828 | 119785
           |  17.76 |   1.29 |   9.50 |   2.11 |   5.98 |   0.57 |  35.07 |  20.86 |  93.12
           |  19.07 |   1.38 |  10.20 |   2.27 |   6.42 |   0.61 |  37.65 |  22.40 |
           |  97.92 |  99.22 |  71.64 |  98.05 |  94.78 |  74.31 |  94.35 |  99.69 |
  ---------+--------+--------+--------+--------+--------+--------+--------+--------+
         1 |    485 |     13 |   4835 |     54 |    423 |    253 |   2700 |     83 |   8846
           |   0.38 |   0.01 |   3.76 |   0.04 |   0.33 |   0.20 |   2.10 |   0.06 |   6.88
           |   5.48 |   0.15 |  54.66 |   0.61 |   4.78 |   2.86 |  30.52 |   0.94 |
           |   2.08 |   0.78 |  28.36 |   1.95 |   5.22 |  25.69 |   5.65 |   0.31 |
  ---------+--------+--------+--------+--------+--------+--------+--------+--------+
  Total       23329     1668    17050     2773     8110      985    47805    26911   128631
              18.14     1.30    13.25     2.16     6.30     0.77    37.16    20.92   100.00


485 PB with 0% of unique features detected
13 PB with 0-25% uniq detected
54 PB with 25-50% uniq detected
423 PB with 50-75% uniq detected
253 PB with 75-100% uniq detected
4835 PB with 100% uniq detected
2700 PB with no uniq detected
83 PB with no expression

*/

/* Make permenant */

data event.refseq_xscripts_with_pacbio;
   set xscript_check;
run;


