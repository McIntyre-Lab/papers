ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Group transcripts based on the proportion of unique features detected into the following:
0%, 0-25%, 25-50%, 50-75%, 75%-100%, 100%

Count the number of transcripts in each bin */

data flag_xscripts_w_unique;
  set event.flag_xscripts_w_unique;
run;

/* Bin transcripts by proportion of unique features detected */

data bin_xscripts;
   set flag_xscripts_w_unique;
   length bin_xscript_perc_uniq_dtct $9.;
   if num_unique_features=0 then bin_xscript_perc_uniq_dtct="no unique";
   else do;
      if perc_unique_features_dtct=. then bin_xscript_perc_uniq_dtct="oops2";
      else if perc_unique_features_dtct=0 then bin_xscript_perc_uniq_dtct="0%";
      else if perc_unique_features_dtct < 0.25 then bin_xscript_perc_uniq_dtct="0-25%";
      else if perc_unique_features_dtct < 0.5 then bin_xscript_perc_uniq_dtct="25-50%";
      else if perc_unique_features_dtct < 0.75 then bin_xscript_perc_uniq_dtct="50-75%";
      else if perc_unique_features_dtct < 1 then bin_xscript_perc_uniq_dtct="75-100%";
      else if perc_unique_features_dtct = 1 then bin_xscript_perc_uniq_dtct="100%";
      else bin_xscript_perc_uniq_dtct="oops";
   end;
run;

/* Check binning */

proc freq data=bin_xscripts;
   tables bin_xscript_perc_uniq_dtct;
run;

/*
  bin_xscript_
  perc_uniq_                               Cumulative    Cumulative
  dtct            Frequency     Percent     Frequency      Percent
  -----------------------------------------------------------------
  0%                 16316       22.19         16316        22.19
  0-25%               1303        1.77         17619        23.96
  100%               12696       17.27         30315        41.23
  25-50%              2045        2.78         32360        44.01
  50-75%              5606        7.62         37966        51.63
  75-100%              669        0.91         38635        52.54
  no unique          34900       47.46         73535       100.00


 Most transcripts have no unique features. For those that do have unique features, most
 have none detected or all detected

*/

/* Make permenant */

data event.xscripts_w_unique_by_bin;
  set bin_xscripts;
run;

