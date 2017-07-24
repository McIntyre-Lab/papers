ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Put all sequence redundancy flags together and decide which events to exclude from further analysis */

data ir_flags;
  set evspl.splicing_events_annot_refseq;
  keep event_id flag_junction_annotated flag_intron_retention;
run;

data chrom_flags;
   set event.rfsq_unannot_hits_flag_chr_gene;
   keep event_id flag_same_chrom flag_same_gene transcript_id;
run;

data event_flags;
   set event.rfsq_unannot_hits_check_event;
   if sum(flag_same_feat1_start,flag_same_feat1_stop,flag_same_feat2_start,flag_same_feat2_stop) > 2;
   keep event_id flag_same_event flag_same_feat1_start flag_same_feat1_stop
                 flag_same_feat2_start flag_same_feat2_stop  transcript_id;
run;

data coord_flags;
   set event.rfsq_unannot_hit_check_coord;
   if sum(flag_feat1_start,flag_feat1_stop,flag_feat2_start,flag_feat2_stop) > 2;
   keep event_id flag_feat1_start flag_feat1_stop  flag_feat2_start flag_feat2_stop  transcript_id;
   rename flag_feat1_start=flag_feat1_start_in_exon
          flag_feat1_stop=flag_feat1_stop_in_exon
          flag_feat2_start=flag_feat2_start_in_exon
          flag_feat2_stop=flag_feat2_stop_in_exon;
run;

data fusion_flags;
  set event.rfsq_unannot_hits_check_ir_fus ;
  keep event_id flag_bad_ir_mult;
  rename flag_bad_ir_mult=flag_ir_from_multigene_fusion;
run;

proc sort data=chrom_flags nodup;
  by event_id  transcript_id;
proc sort data=event_flags nodup;
  by event_id  transcript_id;
proc sort data=coord_flags nodup;
  by event_id  transcript_id;
proc sort data=fusion_flags nodup;
  by event_id;
proc sort data=ir_flags nodup;
  by event_id;
run;


data merge_flags;
  merge chrom_flags event_flags coord_flags;
  by event_id  transcript_id;
  run;

proc sort data=merge_flags nodup;
  by event_id ;
run;

data merge_flags2;
  merge merge_flags (in=in1) fusion_flags (in=in2) ir_flags;
  by event_id;
  if in1 or in2;
run;


data merge_flags3;
   set merge_flags2;
   array change _numeric_ ;
     do over change;
     if change=. then change=0;
   end;
run;

data merge_flags4;
  set merge_flags3;  
  if sum(flag_same_feat1_start,flag_same_feat1_stop,
             flag_same_feat2_start,flag_same_feat2_stop) > 2
             then flag_overlapping_event=1; else flag_overlapping_event=0;
  if sum(flag_feat1_start_in_exon,flag_feat1_stop_in_exon,
             flag_feat2_start_in_exon,flag_feat2_stop_in_exon) = 4
             then flag_event_inside_exon=1; else flag_event_inside_exon=0;
  if sum(flag_feat1_start_in_exon,flag_feat1_stop_in_exon,
             flag_feat2_start_in_exon,flag_feat2_stop_in_exon) ge 3
             then flag_event_overlaps_exon=1; else flag_event_overlaps_exon=0;
run;



/* 
   IR    Same chrom    Same region  same gene   Multigene	OUTCOME
   0     							Ambiguous sequence
   1	 	0						Ambiguous sequence
   1	 	1	0					Ambiguous sequence
   1	 	1	1		0			Ambiguous sequence check multigene
   1	 	1	1		1			Ambiguous sequence within gene
   1	 	1	1		0	1		Ambiguous sequence, multigene region
   1	 	1	1		0	0		Ambiguous sequence

*/


proc freq data=merge_flags4 noprint;
   tables flag_intron_retention*flag_same_chrom*flag_same_gene*flag_overlapping_event*
          flag_event_inside_exon*flag_event_overlaps_exon*flag_ir_from_multigene_fusion / out=flag_summary;
run;

proc print data=flag_summary;
run;

/*

              flag_ flag_                               flag_event_ flag_ir_from_
 flag_intron_ same_ same_ flag_overlapping_ flag_event_  overlaps_    multigene_
   retention  chrom  gene       event       inside_exon     exon        fusion    COUNT

       0        0     0           0              0           0            0          8
       0        1     0           0              0           0            0          2
       0        1     1           0              0           0            0         79
       0        1     1           0              1           1            0          4
       0        1     1           1              0           0            0          9

       1        0     0           0              0           0            0        405
       1        0     0           0              0           0            1          5
       1        1     0           0              0           0            0        205
       1        1     0           0              0           0            1        220
       1        1     0           0              0           1            1         26
       1        1     0           0              1           1            1        454
       1        1     1           0              0           0            0         36
       1        1     1           0              0           0            1         15
       1        1     1           0              0           1            0         35
       1        1     1           0              0           1            1          8


*/

data flag_redundant_seq;
   set merge_flags4;
   length redundancy_class $100.;
   if flag_same_chrom=0 then do;
        if flag_intron_retention=1 then redundancy_class="IR from redundant sequence between chromosomes";
        else redundancy_class="Junction from redundant sequence between chromosomes";
        end;
   else do;
      if flag_same_gene=0 then do;
          if flag_ir_from_multigene_fusion=1 then do;
               if flag_intron_retention=1 then redundancy_class="IR from multigene region";
               else redundancy_class="Junction from multigene region";
               end;
          else do;
               if flag_intron_retention=1 then redundancy_class="IR from different region, same chromosome";
               else redundancy_class="Junction from different region, same chromosome";
               end;
          end;
      else do;
          if flag_intron_retention=1 then redundancy_class="IR from redundant sequence within gene";
          else redundancy_class="Junction from redundant sequence within gene";
      end;
   end;
run;

proc freq data=flag_redundant_seq;
    tables redundancy_class;
run;


/*
 redundancy_class                                      Frequency
 -----------------------------------------------------------------
 IR from different region, same chromosome                  205
 IR from multigene region                                   700
 IR from redundant sequence within gene                      94
 IR from redundant sequence between chromosomes             410
 Junction from different region, same chromosome              2
 Junction from redundant sequence within gene                92
 Junction from redundant sequence between chromosomes         8

*/

/* Make permenant */

data event.refseq_blast_redundancy_class;
   set flag_redundant_seq;
run;




