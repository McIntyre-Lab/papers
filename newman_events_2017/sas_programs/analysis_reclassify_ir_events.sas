ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Reclassifying IR events.

I am going to try the following criteria:
(1) If IR >~ donor exon (90%), then novel donor
(2) Else if IR > ~10 donor then IR
(3) Else unprocessed

Do for IR and for intron and check the concordance */

/* Get intron-to-IR and intron-to-fusion/fragment and merge */
data intron_to_event;
  set event.ir_events_w_intron;
   where flag_no_computed_intron=0; *only keep those with computed introns;
run;

data intron_to_fus_frag;
   set mm10.mm10_introns_to_fragment_fusion;
run;

data ir_coord;
   set evspl.splicing_events_annot_refseq;
   where flag_intron_retention=1;
   keep event_id chr feature1_start feature1_stop
        feature2_start feature2_stop;
run;

proc sort data=intron_to_event;
   by intron_id;
proc sort data=intron_to_fus_frag;
   by intron_id;
run;

data intron_to_ir_fus_frag;
  merge intron_to_event (in=in1) intron_to_fus_frag (in=in2);
  by intron_id;
  if in1 and in2;
run;

proc sort data=intron_to_ir_fus_frag;
  by event_id;
proc sort data=ir_coord;
  by event_id;
run;

data intron_to_ir_W_coord;
  merge intron_to_ir_fus_frag (in=in1) ir_coord (in=in2);
  by event_id;
  if in1 and in2;
  if intron_start = feature2_start-1 then flag_ir_from_3prime=1; else flag_ir_from_3prime=0;
  if intron_stop = feature1_stop-1 then flag_ir_from_5prime=1; else flag_ir_from_5prime=0;
  keep event_id flag_ir_from_3prime flag_ir_from_5prime;
run;

proc freq data=intron_to_ir_w_coord;
   tables flag_ir_from_5prime*flag_ir_from_3prime;
run;

/* First, I want keep only the corresponding donor fusion/fragment for each IR event */

data intron2ir_info;
   set event.intron2ir_frag_fus_all_info;
run;

proc sort data=intron2ir_info;
  by event_id;
proc sort data=intron_to_ir_w_coord;
  by event_id;
run;

data intron2ir_info2;
  merge intron2ir_info (in=in1) intron_to_ir_w_coord (in=in2);
  by event_id;
  if in1 and in2;
run;

data intron2ir_min_info;
   set intron2ir_info2;
   length fusion_id $10.;
   length fragment_id $20.;
   if flag_ir_from_5prime=1 then do;
     fusion_id=fusion_id_5prime;
     fragment_id=fragment_id_5prime;
     flag_fusion_on=flag_fusion_on_5prime;
     mean_apn_fusion=mean_apn_fusion_5prime;
     flag_fragment_on=flag_fragment_on_5prime;
     mean_apn_fragment=mean_apn_fragment_5prime;
     ir_to_fus_ratio=ir_to_5fus_ratio;
     ir_to_frag_ratio=ir_to_5frag_ratio; end;
   else do;
     fusion_id=fusion_id_3prime;
     fragment_id=fragment_id_3prime;
     flag_fusion_on=flag_fusion_on_3prime;
     mean_apn_fusion=mean_apn_fusion_3prime;
     flag_fragment_on=flag_fragment_on_3prime;
     mean_apn_fragment=mean_apn_fragment_3prime;
     ir_to_fus_ratio=ir_to_3fus_ratio;
     ir_to_frag_ratio=ir_to_3frag_ratio; end;
   keep event_id intron_id gene_id flag_no_computed_intron fusion_id fragment_id flag_intron_nsc_on
        mean_apn_intron flag_splicing_on mean_apn_ir flag_fusion_on mean_apn_fusion flag_fragment_on
        mean_apn_fragment ir_to_fus_ratio ir_to_frag_ratio;
run;



data classify_introns_ir;
  set intron2ir_min_info;
   length ir_class_apn1 $20.;
   length ir_class_apn5 $20.;
   length ir_class_apn10 $20.;
   if flag_fusion_on=1 and flag_splicing_on=1 then do; *only do for "analyzable" events;
       /* APN > 1*/
       if mean_apn_ir ge (0.9 * mean_apn_fusion) and mean_apn_ir ge 1 then ir_class_apn1="novel donor";
       else if mean_apn_ir ge (0.1 * mean_apn_fusion) and mean_apn_ir ge 1 then do;
           if mean_apn_intron ge 1 then ir_class_apn1="Probable IR";
           else ir_class_apn1="unprocessed RNA";
           end;
       else ir_class_apn1="unprocessed RNA";

       /* APN > 5*/
       if mean_apn_ir ge (0.9 * mean_apn_fusion) and mean_apn_ir ge 5 then ir_class_apn5="novel donor";
       else if mean_apn_ir ge (0.1 * mean_apn_fusion) and mean_apn_ir ge 5 then do;
           if mean_apn_intron ge 5 then ir_class_apn5="Probable IR";
           else ir_class_apn5="unprocessed RNA";
           end;
       else ir_class_apn5="unprocessed RNA";

       /* APN > 10*/
       if mean_apn_ir ge (0.9 * mean_apn_fusion) and mean_apn_ir ge 10 then ir_class_apn10="novel donor";
       else if mean_apn_ir ge (0.1 * mean_apn_fusion) and mean_apn_ir ge 10 then do;
           if mean_apn_intron ge 10 then ir_class_apn10="Probable IR";
           else ir_class_apn10="unprocessed RNA";
           end;
       else ir_class_apn10="unprocessed RNA";
       end;
    else do;
      ir_class_apn1="low expression";
      ir_class_apn5="low expression";
      ir_class_apn10="low expression";
      end;
run;

proc freq data=classify_introns_ir;
   tables ir_class_apn1 ir_class_apn5 ir_class_apn10;
run;

proc freq data=classify_introns_ir;
   tables flag_splicing_on*flag_fusion_on;
run;

*2745 events where the donor fusion is off but the IR is on;
*34780 events where the donor fusion and IR are on;

/*
 ir_class_apn1      Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
 Probable IR            1502        0.66          1502         0.66
 low expression       193187       84.74        194689        85.40
 novel donor            1180        0.52        195869        85.92
 unprocessed RNA       32098       14.08        227967       100.00


                                             Cumulative    Cumulative
 ir_class_apn5      Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
 Probable IR             265        0.12           265         0.12
 low expression       193187       84.74        193452        84.86
 novel donor             440        0.19        193892        85.05
 unprocessed RNA       34075       14.95        227967       100.00


                                             Cumulative    Cumulative
 ir_class_apn10     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
 Probable IR             103        0.05           103         0.05
 low expression       193187       84.74        193290        84.79
 novel donor             312        0.14        193602        84.93
 unprocessed RNA       34365       15.07        227967       100.00


*/


/* Export APNs and flags so I can make plots of IR vs donor, IR vs intron, intron vs donor */

data data_for_plots;
  set classify_introns_ir;
  keep event_id flag_no_computed_intron mean_apn_intron flag_splicing_on mean_apn_ir
       flag_fusion_on mean_apn_fusion ir_class_apn1 ir_class_apn5 ir_class_apn10;
run;

proc export data=data_for_plots
    outfile="!MCLAB/event_analysis/analysis_output/event_analysis_reclassified_ir.csv"
   dbms=csv replace;
run;


/* Make permenant */

data event.ir_reclassification;
   set classify_introns_ir;
run;
