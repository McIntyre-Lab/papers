ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* RefSeq hits for unannotated junctions : 
Need to figure out if hits represent redundant sequences, or a bug in the creation code */


/* Does the event overlap with an event from the same transcript? */

data flag_gene;
   set event.rfsq_unannot_hits_flag_chr_gene;
run;

data events_by_xs;
  length transcript_id2 $15.;
  set evspl.splicing_events_annot_refseq;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      output;
      end;
  keep event_id transcript_id2 chr strand feature1_start feature1_stop feature2_start feature2_stop;
  rename event_id=transcript_event_id
         transcript_id2=transcript_id
         chr=transcript_chr
         strand=transcript_strand
         feature1_start=transcript_feature1_start
         feature1_stop=transcript_feature1_stop
         feature2_start=transcript_feature2_start
         feature2_stop=transcript_feature2_stop
         ;
run;

proc sort data=events_by_xs;
   by transcript_id;
proc sort data=flag_gene;
   by transcript_id;
run;

%macro mergeXS();

%local iterEvent countEvent;

data _null_;
   set flag_gene nobs=n;
   call symputx('countEvent',n);
   stop;
run;

%let iterEvent=1;

%do %while (&iterEvent. <= &countEvent.);

data event;
   set flag_gene (firstobs=&iterEvent. obs=&iterEvent.);
run;

proc sort data=event;
  by transcript_id;
proc sort data=events_by_xs;
  by transcript_id;
run;

data events_w_junc;
  merge event (in=in1) events_by_xs (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc append base=event_w_xs_junc data=events_w_junc;
run;

%let iterEvent=%eval(&iterEvent.+1);
%end;
%mend;

%mergeXS();

data flag_overlap;
   set event_w_xs_junc;
   if event_id=transcript_event_id then flag_same_event=1;
   else flag_same_event=0; *all should be 0!!;
   if feature1_start=transcript_feature1_start then flag_same_feat1_start=1; else flag_same_feat1_start=0;
   if feature1_stop=transcript_feature1_stop then flag_same_feat1_stop=1; else flag_same_feat1_stop=0;
   if feature2_start=transcript_feature2_start then flag_same_feat2_start=1; else flag_same_feat2_start=0;
   if feature2_stop=transcript_feature2_stop then flag_same_feat2_stop=1; else flag_same_feat2_stop=0;
run;

proc freq data=flag_overlap noprint;
   tables flag_same_event*flag_same_feat1_start*flag_same_feat1_stop
                         *flag_same_feat2_start*flag_same_feat2_stop / out=coord_check;
run;

proc print data=coord_check;
run;

/*
   flag_    flag_same_                  flag_same_
   same_      feat1_      flag_same_      feat2_      flag_same_
   event       start      feat1_stop       start      feat2_stop    COUNT    PERCENT

     0           0             0             0             0         9759    99.2474
     0           0             0             0             1            2     0.0203
     0           0             0             1             1           23     0.2339
     0           0             1             1             1            5     0.0508
     0           1             1             0             0           40     0.4068
     0           1             1             1             0            4     0.0407

As predicted, no events are the same
There are 9 events that could be similar events, perhaps from shorter sequences
However, majority are separate events

*/

/* Make permenant */

data event.rfsq_unannot_hits_check_event;
   set flag_overlap;
run;

