libname event "!MCLAB/event_analysis/sas_data";


/* Identifying annotated and NIC junctions that are detected with STAR but not
   by mapping directly to sequences (event analysis), where detected is defined as
   at least one uniquely-mapping read in at least one replicate.

   I will be checking the following:
    (1) Junction sequence size (is it too small for us to be able to detect?)
    (2) Is the junction sequence not unique?
    (3) How many of these are "annotated" as defined by STAR?
    (3) Is there something about the alignment? (Soft-clipping, spans multiple junctions, etc.)
*/

data check;
   set event.event2star2pacbio_junc_table;
   where flag_star_detected=1 and flag_events_detected=0
   and junction_type ^? "NNC";
   keep junction_id NSC1_apn_star NSC2_apn_star junction_type;
run;

data junc2seq;
  set evspl.mm10_refseq_junc2uniq_seq;
run;

proc sort data=junc2seq;
   by junction_id;
proc sort data=check nodup;
   by junction_id;
run;

data junc2check_seq;
   merge check (in=in1) junc2seq (in=in2);
   by junction_id;
   if in1 and in2;
   junc_seq=length(seq2);
   if junc_seq < 56 then flag_small_event=1;
   else flag_small_event=0;
run;

proc freq data=junc2check_seq;
  tables flag_small_event;
run;

/*                                              Cumulative    Cumulative
               flag_small_event    Frequency     Percent     Frequency      Percent
               ---------------------------------------------------------------------
                              0         127       85.81           127        85.81
                              1          21       14.19           148       100.00



21 of these are sub-read junctions (junctions smaller than the read size)
So why are these kept?
*/

/* How many of these are "annotated" according to STAR? */

data star_junc;
   length chr $5.;
   set event.npc_star_junctions;
   length junction_id $27.;
   length strand_chr $1.;
   if strand=2 then strand_chr="-";
   else strand_chr="+";
   junction_id=catx(":",chr,intron_start-1,intron_stop,strand_chr);
   keep junction_id flag_junction_annotated NSC1 NSC2;
run;

proc sort data=junc2check_seq;
  by junction_id;
proc sort data=star_junc;
  by junction_id;
run;

data junc2check_star;
  merge junc2check_seq (in=in1) star_junc (in=in2);
  by junction_id;
  if in1 and in2;
run;

proc freq data=junc2check_star;
  tables flag_junction_annotated;
run;

/*


   flag_junction_                             Cumulative    Cumulative
        annotated    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0          85       57.43            85        57.43
                1          63       42.57           148       100.00



63 of these are "annotated" in STAR...
*/


/* Make a design file for extracting SAM entries
   Need junction_id sequence sample_to_check small_Event annotated
   coordinates to extract */

data junc_design;
   set junc2check_star;
   retain junction_id seq2 sample flag_small_event flag_junction_annotated chr start stop;
   length sample $4.;
   length chr $5.;
   if NSC1>0 and NSC2>0 then sample="both";
   else if NSC1>0 and NSC2=0 then sample="NSC1";   
   else if NSC2>0 and NSC1=0 then sample="NSC2";   
   else delete;
   chr=scan(junction_id,1,":");
   start=scan(junction_id,2,":")-50;
   stop=scan(junction_id,2,":")+1;
   intron_length=scan(junction_id,3,":")-scan(junction_id,2,":");
   keep junction_id seq2 sample flag_small_event flag_junction_annotated chr start stop intron_length unique_junc_id;
run;

proc sort data=junc_design nodup;
  by _all_;
run;

/* Make FASTA seq for exporting */

data junc_fasta;
  set junc_design;
  length seq_name $81.;
  junc_num=_n_;
  seq_name=catt(">",unique_junc_id);
run;

data junc;
  set junc_fasta;
  keep seq_name junc_num;
  rename seq_name=fasta_seq;
run;

data seq;
  set junc_fasta;
  keep seq2 junc_num;
  rename seq2=fasta_seq;
run;

data junc_fasta2export;
  set junc (in=in1) seq (in=in2);
  if in1 then order=1;
  if in2 then order=2;
run;

proc sort data=junc_fasta2export nodup;
   by junc_num order;
run;

data junc_fasta2export2;
   set junc_fasta2export;
   keep fasta_seq;
run;

/* Export */
proc export data=junc_fasta2export2 outfile="!MCLAB/event_analysis/analysis_output/star_junctions_to_blast.fa"
     dbms=tab replace;
     putnames=no;
run;

proc export data=junc_design outfile="!MCLAB/event_analysis/analysis_output/star_junctions_to_check.csv"
     dbms=csv replace;
run;

