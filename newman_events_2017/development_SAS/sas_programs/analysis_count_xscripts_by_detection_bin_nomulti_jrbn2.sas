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


     data WORK.SPLICE_MATCH    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/references/pacbio_isoforms_list_for_import.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat pacbio_id $10. ; informat source $6.;    informat feature_type $10. ;
        informat start best32. ;  informat end best32.;   informat frame best32. ;
        informat strand $1. ;     informat score best32.; informat transcript_id $18. ;
        informat match_type $23.; informat note $28. ;
        format pacbio_id $10. ;   format source $6. ;     format feature_type $10. ;
        format start best12. ;    format end best12. ;    format frame best12. ;
        format strand $1. ;       format score best12. ;  format transcript_id $18. ;
        format match_type $23. ;  format note $28. ;
        input pacbio_id $ source $ feature_type $ start end
              frame strand $ score transcript_id $ match_type $ note $
              ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data pb2refseq;
   set splice_match;
   where transcript_id ^? "ENSMUST";
   length splice_match_id $18.;
   splice_match_id=scan(transcript_id,1,".");
   keep pacbio_id splice_match_id match_type;
   rename splice_match_id=transcript_id;
run; *8722;

data pb2keep;
  set event.pacbio2refseq_id_nomulti;
  keep pacbio_id;
run;

proc sort data=pb2keep nodup;
  by pacbio_id;
proc sort data=pb2refseq nodup;
  by pacbio_id transcript_id;
run;

data xs_w_pb;
  merge pb2refseq (in=in1) pb2keep (in=in2);
  by pacbio_id;
  if in1 and in2;
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
        0 |      4 |   1307 |     11 |     26 |    408 |   1756
          |   0.07 |  22.23 |   0.19 |   0.44 |   6.94 |  29.86
          |   0.23 |  74.43 |   0.63 |   1.48 |  23.23 |
          |  16.67 |  28.80 |  42.31 |  18.06 |  35.54 |
 ---------+--------+--------+--------+--------+--------+
        1 |     20 |   3231 |     15 |    118 |    740 |   4124
          |   0.34 |  54.95 |   0.26 |   2.01 |  12.59 |  70.14
          |   0.48 |  78.35 |   0.36 |   2.86 |  17.94 |
          |  83.33 |  71.20 |  57.69 |  81.94 |  64.46 |
 ---------+--------+--------+--------+--------+--------+
 Total          24     4538       26      144     1148     5880
              0.41    77.18     0.44     2.45    19.52   100.00
*/

