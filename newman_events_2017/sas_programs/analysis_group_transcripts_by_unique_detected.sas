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
 0%                 23329       22.93         23329        22.93
 0-25%               1668        1.64         24997        24.57
 100%               17050       16.76         42047        41.34
 25-50%              2773        2.73         44820        44.06
 50-75%              8110        7.97         52930        52.03
 75-100%              985        0.97         53915        53.00
 no unique          47805       47.00        101720       100.00

 Most transcripts have no unique features. For those that do have unique features, most
 have none detected or all detected

*/

/* Make permenant */

data event.xscripts_w_unique_by_bin;
  set bin_xscripts;
run;

