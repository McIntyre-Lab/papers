ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/*
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
  set event.rfsq_unannot_hit_check_coord;
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
   set event.rfsq_unannot_hit_check_fusion;
   where flag_multigene=1 and flag_intron_retention=1;
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


/*
 flag_bad_ir_mult
           flag_event2xs_mult

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         1 |     14 |    169 |    183
           |   7.65 |  92.35 | 100.00
           |   7.65 |  92.35 |
           | 100.00 | 100.00 |
  ---------+--------+--------+
  Total          14      169      183
               7.65    92.35   100.00

*/

/* Make permentant */

data event.rfsq_unannot_hits_check_ir_fus;
   set bad_ir_check;
run;

