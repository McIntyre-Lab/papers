ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* RefSeq hits for unannotated junctions : 
Need to figure out if hits represent redundant sequences, or a bug in the creation code */


/*
Check: if on the same chromosome, check if the coordinates of the event overlap with the exonic coordinates of
   the transcript.
*/

data flag_gene;
   set event.rfsq_unannot_hits_flag_chr_gene;
run;

data xs_exons;
   length transcript_id2 $15.;
   set mm10.mm10_exons_w_info;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
        transcript_id2=scan(transcript_id,i,"|");
        output; end;
   keep chrom start stop strand exon_id transcript_id2;
   rename transcript_id2=transcript_id
          chrom=exon_chrom
          start=exon_start
          stop=exon_stop
          strand=exon_strand;
run;

proc sort data=xs_exons nodup;
  by transcript_id;
run;

%macro mergeExons();

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
proc sort data=xs_exons;
  by transcript_id;
run;

data events_w_exons;
  merge event (in=in1) xs_exons (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc append base=event_w_xs_exons data=events_w_exons;
run;

%let iterEvent=%eval(&iterEvent.+1);
%end;
%mend;

%mergeExons();

/* check if the coordinates of the event overlap with the exonic coordinates of the transcript. */

data flag_coords;
   set event_w_xs_exons;
   if flag_same_chrom=1 then do;
      if feature1_start ge exon_start and feature1_start le exon_stop then flag_feat1_start=1; else flag_feat1_start=0;
      if feature1_stop ge exon_start and feature1_stop le exon_stop then flag_feat1_stop=1; else flag_feat1_stop=0;
      if feature2_start ge exon_start and feature2_start le exon_stop then flag_feat2_start=1; else flag_feat2_start=0;
      if feature2_stop ge exon_start and feature2_stop le exon_stop then flag_feat2_stop=1; else flag_feat2_stop=0;
      end;
run;

proc freq data=flag_coords noprint;
   where flag_same_chrom=1;
   tables flag_feat1_start*flag_feat1_stop*flag_feat2_start*flag_feat2_stop*flag_intron_retention*flag_junction_annotated / out=coord_check;
run;

proc print data=coord_check;
run;


/*
 flag_     flag_     flag_     flag_
feat1_    feat1_    feat2_    feat2_    flag_intron_    flag_junction_
 start     stop      start     stop       retention        annotated      COUNT

   0         0         0         0            0                0            215
   0         0         0         0            1                0           5336
   0         0         0         1            0                0              2
   0         0         0         1            1                0              9
   0         0         1         0            0                0              4
   0         0         1         1            0                0             57
   0         1         0         0            0                0              5
   0         1         1         1            1                0             55
   1         0         0         0            0                0              3
   1         0         0         0            1                0              6
   1         1         0         0            0                0             47
   1         1         0         0            1                0             29
   1         1         1         0            1                0             14
   1         1         1         1            0                0              4
   1         1         1         1            1                0            454


*/

/* might not be useful, but collapse event*transcript feature flags */

proc sort data=flag_coords;
   by event_id transcript_id exon_id;
proc means data=flag_coords noprint;
   by event_id transcript_id;
   var flag_feat1_start flag_feat1_stop flag_feat2_start flag_feat2_stop flag_intron_retention;
   output out=flag_coords_by_xs max=;
run;


proc freq data=flag_coords_by_xs noprint;
   tables flag_feat1_start*flag_feat1_stop*flag_feat2_start*flag_feat2_stop*flag_intron_retention / out=coord_check2;
run;

proc print data=coord_check2;
run;

/*

 flag_     flag_     flag_     flag_
feat1_    feat1_    feat2_    feat2_    flag_intron_
 start     stop      start     stop       retention     COUNT

   .         .         .         .            0            8
   .         .         .         .            1          410
   0         0         0         0            0           10
   0         0         0         0            1          436
   0         0         0         1            0            1
   0         0         0         1            1            8
   0         0         1         1            0           24
   0         1         1         1            0            5   ** possible problem
   0         1         1         1            1           52   ** possible problem
   1         0         0         0            0            1
   1         0         0         0            1            2
   1         0         0         1            1            1
   1         0         1         1            0            2   ** possible problem
   1         1         0         0            0           16
   1         1         0         0            1           29
   1         1         0         1            0            1   ** possible problem
   1         1         1         0            0            4   ** possible problem
   1         1         1         0            1           14   ** possible problem
   1         1         1         1            0           30 **** problem here!!
   1         1         1         1            1          457 **** problem here!!



Remainder are probably redundant sequences

So, I estimate that 52+5+2+1+14+14+487=565 are likely due to overlapping coordinates

Of these are there are: 
52+14+457=523 possibly problematic IR events
5+2+1+4+30=42 possible problematic unannotated junctions

However, the unannotated junctions are likely due to sequence ambiguity (since they are split)
*/


/* One last thing to check */

data flag_coords2;
  set flag_coords;
  where sum(flag_feat1_start,flag_feat1_stop,flag_feat2_start,flag_feat2_stop) ge 3 and flag_same_chrom=1;
run;

proc sort data=flag_coords2;
   by gene_id event_id transcript_id feature1_start feature1_stop feature2_start feature2_stop exon_id;
proc means data=flag_coords2 noprint;
   by gene_id event_id transcript_id feature1_start feature1_stop feature2_start feature2_stop ;
   var flag_feat1_start flag_feat1_stop flag_feat2_start flag_feat2_stop flag_intron_retention;
   output out=flag_coords_by_xs2 max=;
run;


proc freq data=flag_coords_by_xs2 noprint;
   tables flag_feat1_start*flag_feat1_stop*flag_feat2_start*flag_feat2_stop*flag_intron_retention / out=coord_check2;
run;

proc print data=coord_check2;
run;

/*
  flag_     flag_     flag_     flag_
 feat1_    feat1_    feat2_    feat2_    flag_intron_
  start     stop      start     stop       retention     COUNT    PERCENT

    0         1         1         1            1           55     10.4364
    1         1         1         0            1           14      2.6565
    1         1         1         1            0            4      0.7590
    1         1         1         1            1          454     86.1480

527 problematic events. Junctions are possibly okay, but the IRs are a problem
*/

/* I want to now check to see how many of these IR events are completely inside a single fusion that contain exons
   from both the event gene and the transcript gene

   I suspect most of these are multi-gene conflicts.
*/

/* Make permenant */

data event.rfsq_unannot_hit_check_coord;
   set flag_coords;
run;

