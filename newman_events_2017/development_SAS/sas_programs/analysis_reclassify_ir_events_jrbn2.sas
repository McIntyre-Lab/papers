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
   set event.mm10_introns_to_fragment_fusion;
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


/*

(1) If IR APN > 5 and IR  ~ exon, then flag_possible_novel_donor=1
(2) If IR APN > 5, IR ~ intron, then flag_possible_IR
(3) If IR APN > 5, IR > intron, then flag_unprocessed_rna=1
(4) If IR APN < 5 then flag_low_expression=1

Figure 1. Possible classifications for putative IR events. Putative IR events are considered representative of a possible novel donor site if the abundance of the event is equivalent to the abundance of its adjacent 5’ exonic region and there is little expression of the adjacent intron. Events are considered as possible IR events if the abundance of the putative IR event and the adjacent intron are equivalent. The remaining set of putative IR events consist of instances where there is low abundance of the event and the intron relative to the 5’ exon, and are considered possible unprocessed transcript.

*/


data classify_introns_ir2;
  set intron2ir_min_info;
   if flag_fusion_on=1 and flag_splicing_on=1 then do; *only do for "analyzable" events;
   flag_low_expressed=0;
       /* APN > 1*/
     if mean_apn_ir ge 5 then do;
       flag_possible_unprocessed=0;
       if mean_apn_ir ge (0.9 * mean_apn_fusion) then flag_possible_novel_donor=1;
       else flag_possible_novel_donor=0;
       if mean_apn_ir ge (0.1 * mean_apn_fusion) then do;
              if mean_apn_intron ge (0.5 * mean_apn_IR) then do;
                    flag_possible_ir=1;
                    flag_possible_unprocessed=0; end;
              else do;
                    flag_possible_ir=0; 
                    flag_possible_unprocessed=1; end;
        end;
        else do;
          flag_possible_ir=0;
          flag_possible_novel_donor=0;
          flag_possible_unprocessed=1; end;
         end;
     else  do;
          flag_possible_ir=0;
          flag_possible_novel_donor=0;flag_possible_unprocessed=1;
     end;  end;
   else flag_low_expressed=1;
run;

proc freq data=classify_introns_ir2 noprint;
   where flag_fusion_on=1 and flag_splicing_on=1;
   tables flag_low_expressed*flag_possible_novel_donor*flag_possible_ir*flag_possible_unprocessed / out=ir_class;
run;

proc print data=ir_class;
run;

data check;
   set classify_introns_ir2;
   where flag_possible_novel_donor=0 and flag_possible_ir=0 and flag_possible_unprocessed=0;
run;


/*
                                    flag_
   flag_low_    flag_possible_    possible_    flag_possible_
   expressed      novel_donor         ir         unprocessed     COUNT    PERCENT

       0               0              0               1          24794    97.7990
       0               0              1               0            200     0.7889
       0               1              0               1            314     1.2386
       0               1              1               0             44     0.1736


       
*/



proc freq data=classify_introns_ir2;
   tables flag_splicing_on*flag_fusion_on;
run;

*1669 events where the donor fusion is off but the IR is on;
*25352 events where the donor fusion and IR are on;


/* Export APNs and flags so I can make plots of IR vs donor, IR vs intron, intron vs donor */

data data_for_plots;
  set classify_introns_ir2;
  keep event_id flag_no_computed_intron mean_apn_intron flag_splicing_on mean_apn_ir
       flag_fusion_on mean_apn_fusion flag_low_expressed flag_possible_novel_donor flag_possible_ir flag_possible_unprocessed ;
run;

proc export data=data_for_plots
    outfile="!MCLAB/event_analysis/analysis_output/event_analysis_reclassified_ir_jrbn3.csv"
   dbms=csv replace;
run;


/* Make permenant */

data event.ir_reclassification_v3;
   set classify_introns_ir2;
run;
