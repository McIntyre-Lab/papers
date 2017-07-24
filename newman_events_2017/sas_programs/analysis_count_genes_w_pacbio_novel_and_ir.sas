ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count genes with different IR classes and unannotated junctions and merge with pacbio data */


data ir_class;
   set event.ir_reclassification_v2;
   where flag_low_expressed = 0;
   length ir_class $50.;
   if flag_possible_novel_donor=1 and flag_possible_ir=1 then ir_class="ambiguous IR (novel donor, IR)";
   else if flag_possible_novel_donor=1 and flag_possible_ir=0 then ir_class="possible novel donor";
   else if flag_possible_novel_donor=0 and flag_possible_ir=1 then ir_class="possible IR";
   else ir_class="possible unprocessed transcript";
   keep event_id ir_class;
run;

proc freq data=ir_class;
  tables ir_class;
run;

data event2gene;
  set evspl.splicing_events_annot_refseq;
  where flag_junction_annotated=0;
  keep event_id gene_id;
run;

proc sort data=event2gene;
  by event_id;
proc sort data=ir_class;
  by event_id;
run;

data ir_class_w_gene;
  merge ir_class (in=in1) event2gene (in=in2);
  by event_id;
  if in1 and in2;
run;

proc sort data=ir_class_w_gene;
   by gene_id ir_class;
proc freq data=ir_class_w_gene noprint;
   by gene_id;
   tables ir_class / out = ir_by_gene;
run;


proc transpose data=ir_by_gene out=ir_class_sbys;
   by gene_id;
   id ir_class;
   var count;
run;


data ir_class_sbys2;
   set ir_class_sbys;
   array change _numeric_;
            do over change;
            if change=. then change=0;
            end;
   drop _NAME_ _LABEL_;
run;

/* Merge with PacBio data */

data pb2refseq;
   set event.pacbio_gene_to_refseq;
run;

proc sort data=pb2refseq;
  by gene_id;
proc sort data=ir_class_sbys2;
  by gene_id;
run;

data pb2ir;
  merge ir_class_sbys2 (in=in1) pb2refseq (in=in2);
  by gene_id;
  if in1 and in2 then flag_has_pacbio=1;
  else if in1 then flag_has_pacbio=0;
  else delete;
run;

data flag_ir;
  set pb2ir;
  if possible_unprocessed_transcript > 0 then flag_unproc=1; else flag_unproc=0;
  if possible_novel_donor > 0 then flag_donor=1; else flag_donor=0;
  if possible_IR > 0 then flag_ir=1; else flag_ir=0;
  if ambiguous_IR__novel_donor__ir_ > 0 then flag_ambig_ir=1; else flag_ambig_ir=0;
  if full_splice_match > 0 and full_splice_match ne . then flag_full_splice=1; else flag_full_splice=0;

  if novel_not_in_catalog > 0 and novel_not_in_catalog ne . then flag_novel_no_cat=1; else flag_novel_no_cat=0;

  if novel_in_catalog > 0 and novel_in_catalog ne . then flag_novel_cat=1; else flag_novel_cat=0;

  if incomplete_splice_match > 0 and incomplete_splice_match ne . then flag_incmp_splice=1; else flag_incmp_splice=0;

  if antisense > 0 and antisense ne . then flag_antisense=1; else flag_antisense=0;

  if genic_intron > 0 and genic_intron ne . then flag_genic_intron=1; else flag_genic_intron=0;

  if genic > 0 and genic ne . then flag_genic=1; else flag_genic=0;

  if intergenic > 0 and intergenic ne . then flag_intergenic=1; else flag_intergenic=0;

  if fusion > 0 and fusion ne . then flag_fusion=1; else flag_fusion=0;
run;




proc freq data=flag_ir noprint;
   tables flag_donor*flag_ambig_ir*flag_ir*flag_unproc / out=ir_counts;
run;


/* 
flag_      flag_                 flag_
donor    ambig_ir    flag_ir    unproc    COUNT

  0          0          0          1      10387
  0          0          1          0         13
  0          0          1          1        223
  0          1          0          0         10
  0          1          0          1         39
  0          1          1          1          8
  1          0          0          0        181
  1          0          0          1        184
  1          0          1          1         11
  1          1          0          1          3
*/

proc print data=ir_counts;
run;


proc freq data=flag_ir noprint;
   tables flag_has_pacbio*flag_donor*flag_ambig_ir*
          flag_novel_no_cat*flag_novel_cat / out=ir_counts;
run;


proc print data=ir_counts;
run;

          
/*
                                   flag_     flag_
flag_has_    flag_      flag_     novel_    novel_
  pacbio     donor    ambig_ir    no_cat      cat     COUNT    PERCENT

    0          0          0          0         0       5114    46.2429
    0          0          1          0         0         25     0.2261
    0          1          0          0         0        285     2.5771
    0          1          1          0         0          1     0.0090
    1          0          0          0         0       3267    29.5415
    1          0          0          0         1        910     8.2286
    1          0          0          1         0        736     6.6552
    1          0          0          1         1        596     5.3893
    1          0          1          0         0         13     0.1176
    1          0          1          0         1          5     0.0452
    1          0          1          1         0          4     0.0362
    1          0          1          1         1         10     0.0904
    1          1          0          0         0         52     0.4702
    1          1          0          0         1         10     0.0904
    1          1          0          1         0         13     0.1176
    1          1          0          1         1         16     0.1447
    1          1          1          0         0          1     0.0090
    1          1          1          0         1          1     0.0090

*/




