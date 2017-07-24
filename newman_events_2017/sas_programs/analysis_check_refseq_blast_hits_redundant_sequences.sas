
ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* RefSeq hits for unannotated junctions : 
Need to figure out if hits represent redundant sequences, or a bug in the creation code

If event and target transcript are on different chromosomes, then flag the sequence as redundant
If on the same chromosome, check if the gene of the event and the gene of the target transcript overlap
If so, and the event is an IR event, check if the complete IR coordinate overlaps with an exonic region
If the event is a junction, then check if the event sequence overlap with junctions from the target gene
*/

/* Chromosome check */

data refseq_hits;
   set event.refseq_blast_hits_w_annot;
   where flag_junction_annotated=0 and flag_unannotated_blast_hit=1;
   keep event_id blast_only_transcript_id;
run;

data refseq_hits2;
   length transcript_id $15.;
   set refseq_hits;
    do i=1 by 1 while(scan(blast_only_transcript_id,i,"|") ^= "");
         transcript_id=scan(blast_only_transcript_id,i,"|");
         output;
         end;
   keep event_id transcript_id;
run;

data event_chrom;
  set evspl.splicing_events_annot_refseq;
  if transcript_id ne "" then flag_junction_annotated=1;
  keep event_id gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop
       flag_junction_annotated flag_intron_retention;
run;

data xs_chrom;
   length transcript_id2 $15.;
   set mm10.mm10_exons_w_info;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
   keep transcript_id2 gene_id chrom;
   rename transcript_id2=transcript_id chrom=transcript_chr gene_id=transcript_gene;
run;

proc sort data=event_chrom;
   by event_id;
proc sort data=refseq_hits2 nodup;
   by event_id;
run;

data refseq_hits_w_event_chr;
  merge refseq_hits2 (in=in1) event_chrom (in=in2);
   by event_id ;
   if in1 and in2;
run;

proc sort data=refseq_hits_w_event_chr nodup;
   by transcript_id;
proc sort data=xs_chrom nodup;
   by transcript_id;
run;

data refseq_hits_w_chroms;
   merge refseq_hits_w_event_chr (in=in1) xs_chrom (in=in2);
   by transcript_id;
   if in1 and in2;
run;

data flag_chrom;
   set refseq_hits_w_chroms;
   if chr=transcript_chr then flag_same_chrom=1;
   else flag_same_chrom=0;
run;

proc freq data=flag_chrom;
   tables flag_same_chrom;
run;

/*
                                              Cumulative    Cumulative
  flag_same_chrom    Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------------
                0         418       27.66           418        27.66
                1        1093       72.34          1511       100.00

*/

/* Check #2: are the genes the same? */

data flag_gene;
  set flag_chrom;
  if transcript_gene=gene_id then flag_same_gene=1;
  else flag_same_gene=0;
run;


proc freq data=flag_gene;
   tables flag_same_gene flag_same_gene*flag_same_chrom;
run;

/*
                              The FREQ Procedure

                                                Cumulative    Cumulative
     flag_same_gene    Frequency     Percent     Frequency      Percent
     -------------------------------------------------------------------
                  0        1325       87.69          1325        87.69
                  1         186       12.31          1511       100.00

           flag_same_gene     flag_same_chrom

           Frequency|
           Percent  |
           Row Pct  |
           Col Pct  |       0|       1|  Total
           ---------+--------+--------+
                  0 |    418 |    907 |   1325
                    |  27.66 |  60.03 |  87.69
                    |  31.55 |  68.45 |
                    | 100.00 |  82.98 |
           ---------+--------+--------+
                  1 |      0 |    186 |    186
                    |   0.00 |  12.31 |  12.31
                    |   0.00 | 100.00 |
                    |   0.00 |  17.02 |
           ---------+--------+--------+
           Total         418     1093     1511
                       27.66    72.34   100.00

*/

/* Check #3: does the event overlap with an event from the transcript? */

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

/*
Check #4: if on the same chromosome, check if the coordinates of the event overlap with the exonic coordinates of
   the transcript.
*/

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
d
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

proc sort data=flag_coords

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
                0           2        0.42             2         0.42
                1         471       99.58           473       100.00

Okay most of these are IR events in multigene fusions. The two that aren't are due to directly adjacent
but non-overlapping exons (ie, exon 1 is from position 100-200, exon 2 is positions 201-300).

I will do one last check (which I have technically already done), whereby for the set of erroneous IR events
I check if the donor exon is in a multigene fusion and flag if so

*/

data ir2exon;
  set evspl.splicing_events_annot_refseq;
  length exon_id $250.;
  where flag_intron_retention=1;
  if index(feature1_id,"intron") > 0 then exon_id=feature2_id;
  else exon_id=feature1_id;
  keep event_id exon_id feature1_id feature2_id;
run;

data fus2exon;
  set mm10.mm10_refseq_fusion_si_info_v2;
  keep fusion_id exon_id flag_multigene;
run;

proc sort data=fus2exon nodup;
   by exon_id;
proc sort data=ir2exon;
   by exon_id;
run;

data ir2fus;
  merge ir2exon (in=in1) fus2exon (in=in2);
  by exon_id;
  if in1 and in2;
run;

data bad_ir;
  set flag_coords;
  where sum(flag_feat1_start,flag_feat1_stop,flag_feat2_start,flag_feat2_stop) ge 3 and flag_same_chrom=1 and flag_intron_retention=1;
run;

proc sort data=bad_ir nodup;
   by event_id;
proc sort data=ir2fus;
   by event_id;
run;

data bad_ir2fus;
  merge bad_ir (in=in1) ir2fus (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=bad_ir2fus;
  tables flag_multigene;
run;


/*

                                               Cumulative    Cumulative
    flag_multigene    Frequency     Percent     Frequency      Percent
    -------------------------------------------------------------------
                 0          35        6.69            35         6.69
                 1         488       93.31           523       100.00


FLAGS:

   IR    Same chrom    Same region  same gene   Multigene	OUTCOME
   0     							Ambiguous sequence
   1	 	0						Ambiguous sequence
   1	 	1	0					Ambiguous sequence
   1	 	1	1		0			Ambiguous sequence check multigene
   1	 	1	1		1			Ambiguous sequence within gene
   1	 	1	1		0	1		Ambiguous sequence, multigene region
   1	 	1	1		0	0		Ambiguous sequence

*/

data bad_ir_multigene;
   set bad_ir2fus;
   where flag_multigene=1;
   keep event_id;
run;

data Event2fus_multi;
   set trim_Event2fus2;
   where flag_multigene=1;
   keep event_id;
run;

proc sort data=bad_ir_multigene nodup;
  by event_id;
proc sort data=event2fus_multi  nodup;
  by event_id;
run;

data bad_ir_check;
  merge bad_ir_multigene (in=in1) event2fus_multi (in=in2);
  by event_id;
  if in1 then flag_bad_ir_mult=1; else flag_bad_ir_mult=0;
  if in2 then flag_event2xs_mult=1; else flag_event2xs_mult=0;
run;

proc freq data=bad_ir_check;
  tables flag_bad_ir_mult*flag_event2xs_mult;
run;


/* Need to consolidate all this!! */



