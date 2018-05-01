/* Libraries */
libname event "!MCLAB/event_analysis/sas_data";
libname evspl "!MCLAB/conesa_pacbio/sas_data/splicing";
libname eventloc "/mnt/store/event_sandbox/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";
ods listing; ods html close;

/* Pull out the 20 NIC junctions that are detected in STAR but not by EA
   Need to figure out what is going on with these: why do we have alignment in one
   but not the other?
   Look at: abundance, junction sequence (is it not unique?), alignment parameters/algorithm, etc. 
*/

data junc2check;
   set event.event2star2pacbio_junc_table;
   where flaG_in_pacbio=1 and flag_in_star=1 and flag_in_catalog=1 and junction_type="NIC/Unannotated";
    if  flag_events_Detected=0 and flag_star_detected=1;
run;

/* no coverage in catalog, very low and inconsistent coverage in STAR (looks like only 1 read maps) */

proc print data=junc2check (keep=junction_id seq_name);
run;

/*
  Obs    junction_id                    seq_name

    1    chr11:49815341:49823211:+      junction_2377653
    2    chr11:101555875:101562296:+    junction_2788723
    3    chr12:84308432:84313936:+      junction_2492648
    4    chr13:49653131:49653528:+      junction_1470580
    5    chr15:79034061:79035256:+      junction_473259
    6    chr16:38548581:38550082:-      junction_1348227
    7    chr18:38297205:38298308:+      junction_2016831
    8    chr19:40340750:40349898:-      junction_2400332
    9    chr2:155335423:155345647:-     junction_2524415
   10    chr3:96696531:96699992:+       junction_448905
   11    chr3:90490272:90490831:-       junction_1072802
   12    chr4:107071018:107080890:+     junction_2170595
   13    chr4:154971554:154972468:+     junction_923646
   14    chr6:47519212:47520792:+       junction_1774982
   15    chr6:30508842:30509684:-       junction_2049198
   16    chr6:30508842:30509684:-       junction_2049198
   17    chr6:122588172:122596631:-     junction_2395616
   18    chr7:101139196:101186387:+     junction_880299
   19    chr8:84864552:84867548:+       junction_789581
   20    chr9:72662613:72677314:+       junction_449174

*/

data junc2seq;
   set evspl.mm10_Refseq_junc2uniq_seq;
run;

proc sort data=junc2check;
   by junction_id;
proc sort data=junc2seq;
   by junction_id;
run;

data junc2check_seq;
  merge junc2check (in=in1) junc2seq (in=in2);
  by junction_id;
  if in1 and in2;
  seq_length=length(seq2);
run;

proc print data=junc2check_seq (keep=junction_id seq_name seq2 seq_length);
run;

/* all junctions are larger than the read size, so it's not a length issue.
   Check if events are from genes with multigene regions, or are in genes that
   aren't expressed (check fusion expression)

   If these aren't in multigene regions, then might try BLASTing these sequences against
   the complete junction catalog to double check there aren't duplicated elsewhere */

/* Get genes expressed and genes with multigene fusion (we dropped these) */

data genes2check;
  set junc2check_seq;
  length gene_id $11.;
  gene_id=scan(event_id,1,":");
  keep gene_id junction_id;
run;

data fus2gene;
  set mm10.mm10_refseq_fusion_si_info_v2;
  keep fusion_id primary_gene_id flag_multigene;
  rename primary_gene_id=gene_id;
run;

proc sort data=genes2check nodup;
  by gene_id;
proc sort data=fus2gene nodup;
  by gene_id;
run;

data fus2check;
  merge genes2check (in=in1) fus2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc means data=fus2check noprint;
  by gene_id;
  var flag_multigene;
  output out=multigene_fusions_per_gene(drop=_TYPE_ _FREQ_)
         sum(flag_multigene)=num_multigene_fusions;
run;

proc print data=multigene_fusions_per_gene;
run;

/*

                num_multigene_
 gene_id           fusions

 11810                      0
 14584                      1
 17966                      0
 17999                      0
 20411                      1
 20615                      1
 207278                     0
 224143                     0
 228812                     0
 229615                     1
 230582                     1
 26912                      1
 269614                     0
 26965                      0
 56736                      0
 66590                      1
 70930                      0
 72649                      0
 77219                      0

Only some of these genes have multigene pieces, so this probably does not explain why these were quantified
but not others
*/

data fus_on;
  set event.flag_fusion_on;
  drop flag_fusion_old_on;
run;

proc sort data=fus_on;
  by fusion_id;
proc sort data=fus2check;
  by fusion_id;
run;

data fus2check_on;
  merge fus2check (in=in1) fus_on (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus2check_on;
  by gene_id flag_multigene;
proc means data=fus2check_on noprint;
  by gene_id flag_multigene;
  var flag_fusion_nsc_on;
  output out=fus_on_by_gene sum(flag_fusion_nsc_on)=num_fusions_dtct;
run;

proc print data=fus_on_by_gene (drop=_TYPE_ );
run;

/*
                        flag_              num_fusions_
  gene_id           multigene    _FREQ_        dtct

  11810                     0      10           10
  14584                     0      18           18
  14584                     1       1            1
  17966                     0      26           24
  17999                     0      29           29
  20411                     0      40           40
  20411                     1       1            1
  20615                     0       3            3
  20615                     1       1            1
  207278                    0      21           21
  224143                    0      11           11
  228812                    0      12           12
  229615                    0      15           15
  229615                    1       1            1
  230582                    0      16           16
  230582                    1       1            1
  26912                     0       9            9
  26912                     1       1            1
  269614                    0      19           19
  26965                     0      24           23
  56736                     0      13           13
  66590                     0      12           12
  66590                     1       1            1
  70930                     0      17           17
  72649                     0      16           16
  77219                     0      11           11

Okay, all of these genes are expressed to some degree, so it can't be a spurious alignment issue.

I want to export the sequences of these junctions and BLAST them against the rest of the catalog to check if there 
are instances of similar sequence.

Then, I want to check the SAM files for the reads that are mapping to these junctions with STAR, to see
   (1) what they are and (2) where they are mapping in the catalog (if anywhere)

Then, I want to realign to both the catalog and with STAR allowing for NO mismatches, and see if these junctions
appear again.

*/

/* Export sequences */

data junc_id;
  set junc2check_seq;
  length fasta_line $81.;
  fasta_line=catt(">",unique_junc_id);
  keep fasta_line unique_junc_id;
run;

data junc_seq;
  set junc2check_seq;
  length fasta_line $81.;
  fasta_line=seq2;
  keep fasta_line unique_junc_id;
run;

data junc2check_fasta;
  set junc_id (in=in1) junc_seq (in=in2);
  if in1 then order=1;
  if in2 then order=2;
run;

proc sort data=junc2check_fasta nodup;
  by unique_junc_id order;
run;

data fasta_for_export;
  set junc2check_fasta;
  keep fasta_line;
run;

proc export data=fasta_for_export outfile="!MCLAB/event_analysis/analysis_output/nic_junctions_to_check.fa"
     dbms=tab replace; putnames=no;
run;




