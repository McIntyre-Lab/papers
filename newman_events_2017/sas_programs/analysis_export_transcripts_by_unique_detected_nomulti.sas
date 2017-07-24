ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Merge in data for all transcripts in catalog (ie, after removal of short sequences), transcripts from expressed genes */
   (1) all transcripts (before gene-level filtering)
   (2) all transcripts in expressed genes (should this be the same as above???)
   (3) transcripts after removal of unlikely genes (ie, xscripts with at least 1 unique detected feature)
*/

/* All transcripts in catalog */

data xscripts_catalog;
  set event.uniq_feature_counts_by_xscript;
  if num_unique_features=0 then do;
      flag_xscript_catalog_has_uniq=0;
      perc_xscript_catalog_uniq=0;
      end;
  else do;
      flag_xscript_catalog_has_uniq=1;
      perc_xscript_catalog_uniq=num_unique_features_dtct/num_unique_features;
      end;
  keep transcript_id num_unique_features num_unique_features_dtct
       flag_xscript_catalog_has_uniq perc_xscript_catalog_uniq;
  rename num_unique_features=num_unique_features_catalog
         num_unique_features_dtct=num_unique_features_dtct_catalog;
run;

/* Transcripts retained after expression-filtering */

data xscripts_expressed;
  set event.feature_cnt_by_xscript_gene_exp;
  if num_unique_features=0 then do;
      flag_xscript_exp_has_uniq=0;
      perc_xscript_exp_uniq=0;
      end;
  else do;
      flag_xscript_exp_has_uniq=1;
      perc_xscript_exp_uniq=num_unique_features_dtct/num_unique_features;
      end;
  keep transcript_id num_unique_features num_unique_features_dtct
       flag_xscript_exp_has_uniq perc_xscript_exp_uniq;
  rename num_unique_features=num_unique_features_exp
         num_unique_features_dtct=num_unique_features_dtct_exp;
run;


/* Transcripts after uniqueness filtering */

data xscripts_probable;
  set event.xscripts_w_unique_by_bin;
  * remove transcripts that have unique features, but none are detected;
  * transcripts with no unique features are okay;
  if flag_xscript_has_unique=1 and num_unique_features_dtct=0 then delete;
  keep transcript_id num_unique_features num_unique_features_dtct
       flag_xscript_has_unique perc_unique_features_dtct;
  rename num_unique_features=num_unique_features_prob
         num_unique_features_dtct=num_unique_features_dtct_prob
         flag_xscript_has_unique=flag_xscript_prob_has_uniq
         perc_unique_features_dtct=perc_xscript_prob_uniq;
run;
  

/* Merge together, fill out any empty spaces with zeros (use flags for assigning groups in python) */

proc sort data=xscripts_catalog;
   by transcript_id;
proc sort data=xscripts_expressed;
   by transcript_id;
proc sort data=xscripts_probable;
   by transcript_id;
run;

data xscripts_all;
   merge xscripts_catalog (in=in1) xscripts_expressed (in=in2) xscripts_probable (in=in3);
   by transcript_id;
   if in1 then flag_xscript_in_catalog=1; else flag_xscript_in_catalog=0;
   if in2 then flag_xscript_expressed=1; else flag_xscript_expressed=0;
   if in3 then flag_xscript_probable=1; else flag_xscript_probable=0;
run;

data xscripts_all2;
   set xscripts_all;
   array change _numeric_;
        do over change;
            if change=. then change=0;
        end;
run ;


/* Check to make sure everything is in catalog, and:
   catalog > expressed > probable */

proc freq data=xscripts_all2 noprint;
  tables flag_xscript_in_catalog*flag_xscript_expressed*flag_xscript_probable / out=xscript_cnt;
run;

proc print data=xscript_cnt;
run;

/*

   flag_xscript_    flag_xscript_    flag_xscript_
     in_catalog       expressed         probable      COUNT    PERCENT

         1                0                0          18283    19.9122
         1                1                0          16316    17.7699
         1                1                1          57219    62.3178

*/


/* Export data for making plots to compare between lists of transcripts */

proc export data=xscripts_all2
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_transcripts_for_comparison_nomulti.csv" dbms=csv replace;
run;

/* Export data for making plots to look at transcript binning */

proc export data=event.xscripts_w_unique_by_bin
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_transcripts_by_perc_unique_detected_nomulti.csv" dbms=csv replace;
run;

