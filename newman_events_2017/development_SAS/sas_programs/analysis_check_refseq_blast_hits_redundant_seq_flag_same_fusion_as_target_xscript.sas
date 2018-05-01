ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* RefSeq hits for unannotated junctions : 
Need to figure out if hits represent redundant sequences, or a bug in the creation code */




/* I want to now check to see how many of these IR events are completely inside a single fusion that contain exons
   from both the event gene and the transcript gene

   I suspect most of these are multi-gene conflicts.
*/

data fusions_for_event;
  set mm10.mm10_fusion_2_gene_id;
  keep fusion_id gene_id flag_multigene fusion_start fusion_stop;
  rename fusion_id=event_fusion_id
         fusion_start=event_fusion_start
         fusion_stop=event_fusion_stop;
run;

data fusions_for_xs;
  set mm10.mm10_fusion_2_gene_id;
  keep fusion_id transcript_id fusion_start fusion_stop;
  rename fusion_id=xs_fusion_id
         fusion_start=xs_fusion_start
         fusion_stop=xs_fusion_stop;
run;

proc sort data=fusions_for_event nodup;
   by gene_id;
proc sort data=fusions_for_xs nodup;
   by transcript_id;
run;

proc datasets noprint;
  delete event_w_fusions;
run;

data flag_coords2;
  set  event.rfsq_unannot_hit_check_coord;
  where sum(flag_feat1_start,flag_feat1_stop,flag_feat2_start,flag_feat2_stop) ge 3 and flag_same_chrom=1;
run;


%macro mergeFusion();

%local iterEvent countEvent iterXS countXS;

data _null_;
   set flag_coords2 nobs=n;
   call symputx('countEvent',n);
   stop;
run;

%let iterEvent=1;
%do %while (&iterEvent. <= &countEvent.);

data event;
   set flag_coords2 (firstobs=&iterEvent. obs=&iterEvent.);
run;

proc sort data=event;
  by gene_id;
proc sort data=fusions_for_event;
  by gene_id;
run;

data events_w_fus;
  merge event (in=in1) fusions_for_event (in=in2);
  by gene_id;
  if in1 and in2;
run;

data events_w_fus1;
   set events_w_fus;
   * keep only fusions where the event is contained completely within a fusion;
   if feature1_start ge event_fusion_start and feature2_stop le event_fusion_stop;
run;

proc append base=event_w_fusions data=events_w_fus1;
run;

%let iterEvent=%eval(&iterEvent.+1);

%end;
%mend;

%mergeFusion();

%macro mergeXSFus();
%local  iterXS countXS;
   data _null_ ;
      set event_w_fusions nobs=n;
      call symputx('countXS',n);
      stop;
   run;

    %let iterXS=1;

    %do %while (&iterXS. <= &countXS.);

    data event2;
       set event_w_fusions (firstobs=&iterXS. obs=&countXS.);
    run;

    proc sort data=event2;
      by transcript_id;
    proc sort data=fusions_for_xs;
      by transcript_id;
    run;

    data events_w_fus2;
      merge event2 (in=in1) fusions_for_xs (in=in2);
      by transcript_id;
      if in1 and in2;
    run;

    proc append base=event_w_fusions2 data=events_w_fus2;
    run;

    %let iterXS=%eval(&iterXS.+1);

    %end;
    %mend;

%mergeXSFus();


data trim_event2fus;
  set event_w_fusions2;
  if event_fusion_id=xs_fusion_id;
run;

proc sort data=trim_event2fus nodup out=trim_event2fus2;
   by _all_;
run;

/* Check: are the set of remaining fusions all multigene? */

proc freq data=trim_Event2fus2;
  tables flag_multigene;
run;


/*
                                             Cumulative    Cumulative
  flag_multigene    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------------
               0           2        0.45             2         0.45
               1         441       99.55           443       100.00


Okay most of these are IR events in multigene fusions. The two that aren't are due to directly adjacent
but non-overlapping exons (ie, exon 1 is from position 100-200, exon 2 is positions 201-300). */

/* Make permenant */

data event.rfsq_unannot_hit_check_fusion;
   set trim_Event2fus2;
run;

