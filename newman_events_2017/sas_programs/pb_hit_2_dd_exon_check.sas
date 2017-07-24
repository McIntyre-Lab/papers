ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For the genes that have DD exons (or all events???), how many have hits to PacBio transcripts?

   since I did the BLAST on fragments, I will need to convert fragments to fusions
   (ie, if a fragment has a hit, then so does the fusion)
   ie, exon is DD: ## with PB hits, ## without PB hits

   Also, flag if PB isoform is also DD, and where it is being seen
   Table:

   exon_id, pacbio_id, pb_npc, pb_opc, exon_npc, exon_opc */

data dd_genes;
   set event.mm10_flag_gene_on_by_cell_apn5;
   /* If there are more exons detected in either cell type than in both cell types together,
      then the gene has differentially detected exons */
    if sum((flag_fusion_nsc_on-flag_fusion_both_on),
             (flag_fusion_old_on-flag_fusion_both_on),
             flag_fusion_both_on) = 1 then flag_gene_monoexon=1;
    else flag_gene_monoexon=0;

   if flag_gene_both_on=1 then do;
      flag_gene_cell_specific=0;
      if flag_fusion_nsc_on > flag_fusion_both_on or flag_fusion_old_on > flag_fusion_both_on then flag_gene_exon_dd=1;
       else flag_gene_exon_dd=0;
      end;
   else do;
       if flag_gene_nsc_on=0 and flag_gene_old_on=0 then flag_gene_cell_specific=.;
       else flag_gene_cell_specific=1;
   end;

   if flag_gene_exon_dd=1 then output;
   keep gene_id;
run;

*1560 genes in total;

/* Get the exons of these genes that are detected at APN>5 */

data fus2gene;
  set mm10.mm10_refseq_fusion_si_info_v2;
  keep primary_gene_id fusion_id;
  rename primary_gene_id=gene_id;
run;

data fus_dtct;
  set event.flag_fusion_on_apn5;
run;

proc sort data=fus2gene nodup;
  by gene_id fusion_id;
proc sort data=dd_genes;
  by gene_id;
run;

data fus2gene_dd;
  merge fus2gene (in=in1) dd_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=fus_dtct;
  by fusion_id;
proc sort data=fus2gene_dd;
  by fusion_id;
run;

data fus_dtct_gene_dd;
  merge fus2gene_dd (in=in1) fus_dtct (in=in2);
  by fusion_id;
  if in1 and in2;
run;


/* What pieces have PB hits ? */

data frag2pb;
  set event.blast_pacbio_fragments_in_refseq;
  where flag_fragment_on=1 and flag_frag_redundant=0;
  keep fragment_id pacbio_id;
run;

data frag2fus;
   set mm10.mm10_exon_fragment_flagged;
   keep fragment_id fusion_id;
run;

proc sort data=frag2pb nodup;
  by fragment_id;
proc sort data=frag2fus;
  by fragment_id;
run;

data frag2fus2pb;
  merge frag2fus (in=in1) frag2pb (in=in2);
  by fragment_id;
  if in2 then flag_pb_hit=1; else flag_pb_hit=0;
  if in1 then output;
run;

proc sort data=frag2fus2pb;
  by fusion_id pacbio_id;
proc means data=frag2fus2pb noprint;
  by fusion_id pacbio_id;
  var flag_pb_hit;
  output out=flag_fus_pb_hit max=;
run;

proc sort data=flag_fus_pb_hit;
  by fusion_id;
proc sort data=fus_dtct_gene_dd;
  by fusion_id;
run;

data fus2gene_dd_w_pb;
  merge fus_dtct_gene_dd (in=in1) flag_fus_pb_hit (in=in2);
  by fusion_id;
  if in1 and in2;
  drop _TYPE_ _FREQ_;
run;

/* Flag if fusion is testable for DD, and whether fusion is DD */

data flag_fus_dd;
   set fus2gene_dd_w_pb;
   if flag_fusion_nsc_on=. or flag_fusion_old_on=. then flag_dd_testable=0;
   else if flag_fusion_nsc_on=0 and flag_fusion_old_on=0 then flag_dd_testable=0;
   else do;
     flag_dd_testable=1;
     if flag_fusion_nsc_on=1 and flag_fusion_old_on=1 then flag_fusion_dd=0;
     else flag_fusion_dd=1;
     end;
run;

proc freq data=flag_fus_dd noprint;
  tables flag_dd_testable*flag_fusion_dd*flag_pb_hit / out=dd_check;
run;

proc print data=dd_check;
run;

/* this is actually the fusion*pacbio_id count
               flag_
  flag_dd_    fusion_    flag_pb_
  testable       dd         hit      COUNT

      0          .           0       12972 <-ignore
      0          .           1        2459 <-ignore

      1          0           0        6581
      1          0           1        7967
      1          1           0        3065
      1          1           1        1425

*/

/* Flag if PB transcript is also DD */

* Import lists;

proc import datafile="!MCLAB/event_analysis/text_data/TOTAL_Intersection.txt"
     out=pb_both dbms=tab replace;
     guessingrows=20000;
run;

proc import datafile="!MCLAB/event_analysis/text_data/TOTAL_OPC_exclusive.txt"
     out=pb_opc dbms=tab replace;
     guessingrows=20000;
run;

proc import datafile="!MCLAB/event_analysis/text_data/TOTAL_NPC_exclusive.txt"
     out=pb_npc dbms=tab replace;
     guessingrows=20000;
run;

data pb_both2;
  set pb_both;
  keep transcript_id;
  rename transcript_id=pacbio_id;
run;

data pb_npc2;
  set pb_npc;
  keep transcript_id;
  rename transcript_id=pacbio_id;
run;

data pb_opc2;
  set pb_opc;
  keep transcript_id;
  rename transcript_id=pacbio_id;
run;

proc sort data=pb_both2;
  by pacbio_id;
proc sort data=pb_npc2;
  by pacbio_id;
proc sort data=pb_opc2;
  by pacbio_id;
run;

data pb_as;
  merge pb_both2 (in=in1) pb_npc2 (in=in2) pb_opc2 (in=in3);
  by pacbio_id;
  if in1 then flag_in_both=1; else flag_in_both=0;
  if in2 then flag_in_npc_only=1; else flag_in_npc_only=0;
  if in3 then flag_in_opc_only=1; else flag_in_opc_only=0;
run;

proc freq data=pb_as noprint;
  tables flag_in_both*flag_in_npc_only*flag_in_opc_only / out=pb_count;
proc print data=pb_count;
run;

data pb_as_flag;
  set pb_as;
  if flag_in_both=1 then flag_pacbio_dd=0;
  else flag_pacbio_dd=1;
run;

proc sort data=pb_as_flag;
  by pacbio_id;
proc sort data=flag_fus_dd;
  by pacbio_id;
run;

data pb_vs_fus_dd;
  merge flag_fus_dd (in=in1) pb_as_flag (in=in2);
  by pacbio_id;
  if in1 and in2;
run;

data pb_vs_fus_dd2 ambig;
  set pb_vs_fus_dd;
  if flag_fusion_nsc_on=. or flag_fusion_old_on=. then output ambig;
  else output pb_vs_fus_dd2;
run;

data pb_vs_fus_dd3;
  set pb_vs_fus_dd2;
  if flag_pacbio_dd=0 then flag_pacbio_not_dd=1;
  else if flag_pacbio_dd=1 then flag_pacbio_not_dd=0;
  if flag_fusion_dd=. then delete;
run; *8989 combos to test;

proc sort data=pb_vs_fus_dd3;
   by fusion_id;
proc means data=pb_vs_fus_dd3 noprint;
   by fusion_id;
   var flag_pb_hit flag_fusion_dd flag_pacbio_dd flag_pacbio_not_dd;
   output out=pb2fus_dd max(flag_pb_hit)=flag_fusion_has_pb_hit
                        max(flag_fusion_dd)=flag_fusion_dd
                        max(flag_pacbio_dd)=flag_pacbio_dd
                        sum(flag_pb_hit)=num_pacbio_hits
                        sum(flag_pacbio_dd)=num_pacbio_dd
                        sum(flag_pacbio_not_dd)=num_pacbio_not_dd;
run;
*4177 fusions;

proc freq data=pb2fus_dd noprint;
  tables flag_fusion_has_pb_hit*flag_fusion_dd*flag_pacbio_dd / out=fus_check;
proc print data=fus_check;
run;

/*
                   flag_      flag_
  flag_fusion_    fusion_    pacbio_
   has_pb_hit        dd         dd      COUNT

        1            0          0        2856
        1            0          1         389
        1            1          0         834
        1            1          1          98


*/

data fus2check;
  set pb2fus_dd;
  where flag_fusion_dd=1 and flag_pacbio_dd=0;
  keep fusion_id;
run;

proc sort data=fus2check;
  by fusion_id;
proc sort data=pb_vs_fus_dd3;
  by fusion_id;
run;

data fus2check2;
  merge fus2check (in=in1) pb_vs_fus_dd3 (in=in2);
  by fusion_id;
  if in1 and in2;
run;

/* Import LR counts */

%macro importLR(sample,name);
proc import datafile="!MCLAB/event_analysis/text_data/LR_SR_cor/&sample..5merge.abundance.txt"
     out=&name._lr dbms=tab replace; guessingrows=20000;
run;

data &name._lr2;
  set &name._lr;
  &name._fl_reads=count_fl;
  &name._total_reads=count_fl+count_nfl;
  keep pbid &name._fl_reads &name._total_reads;
  rename pbid=pacbio_id;
run;

proc sort data=&name._lr2;
   by pacbio_id;
run;

%mend;

%importLR(NSC1,npc1);
%importLR(NSC2,npc2);
%importLR(OLD1,opc1);
%importLR(OLD2,opc2);

data all_lr_counts;
  merge npc1_lr2 npc2_lr2 opc1_lr2 opc2_lr2 ;
  by pacbio_id;
run;

data mean_lr_counts;
  set all_lr_counts;
  mean_npc_fl=mean(npc1_fl_reads,npc2_fl_reads);
  mean_npc_all=mean(npc1_total_reads,npc2_total_reads);
  mean_opc_fl=mean(opc1_fl_reads,opc2_fl_reads);
  mean_opc_all=mean(opc1_total_reads,opc2_total_reads);
  diff_fl=mean_opc_fl-mean_npc_fl;
  diff_total=mean_opc_all-mean_npc_all;
run;

data lr_counts;
  set mean_lr_counts;
  keep pacbio_id mean_npc_all mean_opc_all diff_fl diff_total;
run;

proc sort data=lr_counts;
  by pacbio_id;
proc sort data=fus2check2;
  by pacbio_id;
run;

data fus_lr_counts;
  merge fus2check2 (in=in1) lr_counts (in=in2);
  by pacbio_id;
  if in1 and in2;
run;

data flag_diffs;
   set fus_lr_counts;
   if abs(diff_total) > 5 then flag_pacbio_lr_dd=1; else flag_pacbio_lr_dd=0;
   if mean_npc_all>mean_opc_all then flag_pb_higher_in_npc=1; flag_pb_higher_in_npc=0;
   if mean_npc_all<mean_opc_all then flag_pb_higher_in_opc=1; flag_pb_higher_in_opc=0;
run;

proc sort data=flag_diffs;
  by fusion_id;
proc means data=flag_diffs noprint;
  by fusion_id;
  var flag_fusion_nsc_on flag_fusion_opc_on flag_pacbio_lr_dd flag_pb_higher_in_npc flag_pb_higher_in_old;
  output out=fus2pb_dd_check max=;
run;

proc freq data=fus2pb_dd_check noprint;
  tables flag_fusion_nsc_on*flag_fusion_opc_on*flag_pacbio_lr_dd*
         flag_pb_higher_in_npc*flag_pb_higher_in_old / out=fus2pb_count;
run;

proc print data=fus2pb_count;
run;






