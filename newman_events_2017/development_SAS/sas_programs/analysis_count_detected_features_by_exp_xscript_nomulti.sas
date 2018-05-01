ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For probable transcripts, calculate the number of detected features total and by type,
   the number of unique total and by type, and the percent unique total and by type

   I also want to determine the number of transcripts that have unique exons,
   fragments and junctions, and what the overlap is between these sets (for making Venn diagrams) 
   ie. whole exon, exon fragments, splicing events */


data feat2xs2gene;
  set event.feature2xs2gene_exp_only_nomulti;
run;

data feature_data;
   set event.feature_data_for_plots_gene_exp;
   where feature_type ne 'fusion';
run;

proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=feature_data;
   by feature_id;
run;

data feat_data_w_xs;
  merge feat2xs2gene (in=in1) feature_data (in=in2);
  by feature_id;
  if in1 and in2;
run;

proc sort data=feat_data_w_xs;
   by transcript_id;
run;

/* Calculate total detected by type */
* Singletons;
proc means data=feat_data_w_xs noprint;
   where feature_type="fragment" and flag_singleton=1;
   by transcript_id;
   var flag_feature_on;
   output out=num_singleton_per_xs sum=num_singleton_dtct;
run;

* Fragments;
proc means data=feat_data_w_xs noprint;
   where feature_type="fragment" and flag_singleton=0;
   by transcript_id;
   var flag_feature_on;
   output out=num_fragment_per_xs sum=num_fragment_dtct;
run;

* Junctions;
proc means data=feat_data_w_xs noprint;
   where feature_type="splicing";
   by transcript_id;
   var flag_feature_on;
   output out=num_splicing_per_xs sum=num_splicing_dtct;
run;

/* Calculate unique detected by type */
* Singletons;
proc means data=feat_data_w_xs noprint;
   where feature_type="fragment" and flag_singleton=1 and flag_feature_unique=1;
   by transcript_id;
   var flag_feature_on;
   output out=num_uniq_singleton_per_xs sum=num_uniq_singleton_dtct;
run;

* Fragments;
proc means data=feat_data_w_xs noprint;
   where feature_type="fragment" and flag_singleton=0 and flag_feature_unique=1;
   by transcript_id;
   var flag_feature_on;
   output out=num_uniq_fragment_per_xs sum=num_uniq_fragment_dtct;
run;

* Junctions;
proc means data=feat_data_w_xs noprint;
   where feature_type="splicing" and flag_feature_unique=1;
   by transcript_id;
   var flag_feature_on;
   output out=num_uniq_splicing_per_xs sum=num_uniq_splicing_dtct;
run;


/* Merge and calculate perc uniq */
  
proc sort data=num_singleton_per_xs (drop=_TYPE_ _FREQ_);
   by transcript_id;
proc sort data=num_fragment_per_xs (drop=_TYPE_ _FREQ_);
   by transcript_id;
proc sort data=num_splicing_per_xs (drop=_TYPE_ _FREQ_);
   by transcript_id;
proc sort data=num_uniq_singleton_per_xs (drop=_TYPE_ _FREQ_);
   by transcript_id;
proc sort data=num_uniq_fragment_per_xs (drop=_TYPE_ _FREQ_);
   by transcript_id;
proc sort data=num_uniq_splicing_per_xs (drop=_TYPE_ _FREQ_);
   by transcript_id;
run;

data num_features_by_xs;
   merge num_singleton_per_xs num_fragment_per_xs num_splicing_per_xs
         num_uniq_singleton_per_xs num_uniq_fragment_per_xs num_uniq_splicing_per_xs;
   by transcript_id;
run;

data num_features_by_xs2;
   set num_features_by_xs;
   array change _numeric_;
        do over change;
            if change=. then change=0;
        end;
run;

* Calc percent;
data perc_features_uniq;
   set num_features_by_xs2;
   num_total_features_dtct=num_singleton_dtct+num_fragment_dtct+num_splicing_dtct;
   num_uniq_features_dtct=num_uniq_singleton_dtct+num_uniq_fragment_dtct+num_uniq_splicing_dtct;
   /* Calc perc unique for each */
   perc_features_uniq=num_uniq_features_dtct/num_total_features_dtct;
   perc_singleton_uniq=num_uniq_singleton_dtct/num_singleton_dtct;
   perc_fragment_uniq=num_uniq_fragment_dtct/num_fragment_dtct;
   perc_splicing_uniq=num_uniq_splicing_dtct/num_splicing_dtct;
   if num_uniq_singleton_dtct > 0 then flag_xs_has_uniq_singleton=1; else flag_xs_has_uniq_singleton=0;
   if num_uniq_fragment_dtct > 0 then flag_xs_has_uniq_fragment=1; else flag_xs_has_uniq_fragment=0;
   if num_uniq_splicing_dtct > 0 then flag_xs_has_uniq_splicing=1; else flag_xs_has_uniq_splicing=0;
run;

/* Merge in transcript bins, in case I want to split these by list */
data xscript_bins;
   set event.xscripts_w_unique_by_bin;
   keep transcript_id bin_xscript_perc_uniq_dtct;
run;

proc sort data=perc_features_uniq;
  by transcript_id;
proc sort data=xscript_bins;
  by transcript_id;
run;

data perc_features_uniq_w_bin;
  merge xscript_bins (in=in1) perc_features_uniq (in=in2);
  by transcript_id;
  if in1 and in2 then output;
run;


/* Overlap between transcripts resolved by unique piece */
proc freq data=perc_features_uniq_w_bin noprint;
    tables flag_xs_has_uniq_singleton*flag_xs_has_uniq_fragment*flag_xs_has_uniq_splicing
         / out=count_xs_by_reslv_seq;
run;

proc print data=count_xs_by_reslv_seq;
run;

/*
   flag_xs_     flag_xs_     flag_xs_
  has_uniq_    has_uniq_    has_uniq_
  singleton     fragment     splicing    COUNT    PERCENT

      0            0            0        51216    69.6485
      0            0            1         4187     5.6939
      0            1            0         3860     5.2492
      0            1            1         1480     2.0126
      1            0            0         6959     9.4635
      1            0            1         5203     7.0755
      1            1            0          214     0.2910
      1            1            1          416     0.5657



Singleton only	6959
Fragment only	3860
Splicing only	4187
Single+Frag	214
Single+Splice	5203
Frag+Splice	1480
All		416

*/

/* Make permenant and export */

data event.feature_dtct_cnt_by_xscript_exp;
   set perc_features_uniq_w_bin;
run;

proc export data=event.feature_dtct_cnt_by_xscript_exp
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_detected_features_by_type_probable_transcripts_nomulti.csv"
     dbms=csv replace;
run;


