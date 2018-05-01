/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the benchmarked simulation, import the STAR junction results and EA junction coverage counts and check
   to see which method captures the junctions of these genes 

   Note that the catalog-based approach will miss junctions from transcripts we could not identify

   Flag junctions from the identified simulated transcripts

*/

/* First I want to get a list of possible junctions from the transcripts simulated */

data sim_exons_mm10;
   set event.simulated_exons_mm9_to_mm10;
   where flag_mm10_coord=1;
   length sim_gene_id $20.;
   sim_gene_id=scan(sim_exon_id,1,":");
   keep simulation test sim_exon_id sim_gene_id chr exon_start exon_stop strand;
run;

* per gene, count the number of exons. I am going to drop genes that only have a single exon, since these
  won't have junctions;

proc sort data=sim_exons_mm10;
   by simulation test sim_gene_id;
proc freq data=sim_exons_mm10 noprint;
   by simulation test;
   tables sim_gene_id / out=ex_per_gene;
run;

data genes2drop;
  set ex_per_gene;
  if count < 2 then output;
  keep simulation test sim_gene_id;
run;

data sim_exons_mm10_2;
  merge sim_exons_mm10 (in=in1) genes2drop (in=in2);
  by simulation test sim_gene_id;
  if in1 and in2 then delete;
  else if in1 then output;
run;

/* Create annotated junctions for these transcripts/genes */

proc sort data=sim_exons_mm10_2;
   by simulation test sim_gene_id exon_start exon_stop;
run;

data sim_junctions;
   set sim_exons_mm10_2;
   by simulation test sim_gene_id;
   donor_site=lag(exon_stop);
   donor_gene=lag(sim_gene_id); *if this doesn't match with sim_gene_id then delete;
run;

   
data sim_junctions2;
   set sim_junctions;
   length junction_id $50.;
   if donor_gene ^= sim_gene_id then delete;
   junction_id=catx(":",chr,donor_site,exon_start,strand);
   keep simulation test junction_id;
run;

proc sort data=sim_junctions2 nodup;
   by simulation test junction_id;
run;

/* (1) Flag simulated junctions where the RefSeq transcript is identified
   (2) Flag simulated junctions that are in the catalog (and what type of junction they are)
   (3) Merge in APN counts for simulated junctions
   (4) Merge in STAR counts for simulated junctions
   (5) Count overlaps by sample */

/* For each simulated sample, identify the junctions that have an identified RefSeq transcript */

data xscript_s1t1 xscript_s2t1 xscript_s1t2 xscript_s2t2 xscript_s3t1 xscript_s3t2;
   set event.sim_xscript2refseq_blast_best;
   if simulation="sim1" and test="test1" then output xscript_s1t1;
   if simulation="sim2" and test="test1" then output xscript_s2t1;
   if simulation="sim3" and test="test1" then output xscript_s3t1;
   if simulation="sim1" and test="test2" then output xscript_s1t2;
   if simulation="sim2" and test="test2" then output xscript_s2t2;
   if simulation="sim3" and test="test2" then output xscript_s3t2;
   keep transcript_id;
run;

proc sort data=xscript_s1t1 nodup;
  by transcript_id;
proc sort data=xscript_s2t1 nodup;
  by transcript_id;
proc sort data=xscript_s3t1 nodup;
  by transcript_id;
proc sort data=xscript_s1t2 nodup;
  by transcript_id;
proc sort data=xscript_s2t2 nodup;
  by transcript_id;
proc sort data=xscript_s3t2 nodup;
  by transcript_id;
run;

data xscript_sim;
   merge xscript_s1t1 (in=in1) xscript_s2t1 (in=in2) xscript_s3t1 (in=in3) 
         xscript_s1t2 (in=in4) xscript_s2t2 (in=in5) xscript_s3t2 (in=in6);
   by transcript_id;
   if in1 then flag_sim1_test1_xs=1; else flag_sim1_test1_xs=0;
   if in2 then flag_sim2_test1_xs=1; else flag_sim2_test1_xs=0;
   if in3 then flag_sim3_test1_xs=1; else flag_sim3_test1_xs=0;
   if in4 then flag_sim1_test2_xs=1; else flag_sim1_test2_xs=0;
   if in5 then flag_sim2_test2_xs=1; else flag_sim2_test2_xs=0;
   if in6 then flag_sim3_test2_xs=1; else flag_sim3_test2_xs=0;
run;

data event2xs;
   set evspl.splicing_events_annot_refseq;
   length junction_id $50.;
   length transcript_id2 $20.;
   junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
   if transcript_id = "" then delete;
   else do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
        transcript_id2=scan(transcript_id,i,"|");
        output; end;
   keep junction_id transcript_id2;
   rename transcript_id2=transcript_id;
run;

proc sort data=xscript_sim;
  by transcript_id;
proc sort data=event2xs nodup;
  by transcript_id;
run;

data event_sim;
  merge event2xs (in=in1) xscript_sim (in=in2);
  by transcript_id;
  if in1 and in2;
  drop transcript_id;
run;

proc sort data=event_sim;
   by junction_id;
proc means data=event_sim noprint;
   by junction_id;
   var flag_: ;
   output out=event_sim2(drop= _TYPE_ _FREQ_) max=;
run;

/* Merge with the "simulated junctions" list 
   For a junction in a given sample, if the corresponding "flag_*_xs" is 1, that means that
   a transcript that contains that junction was identiable via BLAST for that sample.
   Otherwise, a corresponding transcript was not identified 

   E.g. for sim1_test1, if flag_sim1_test1=1 then flag_transcript_identified=1
   else flag_transcript_identified=0 */


proc sort data=event_sim2;
   by junction_id;
proc sort data=sim_junctions2;
   by junction_id;
run;

data sim_junctions_flagged;
   merge sim_junctions2 (in=in1) event_sim2 (in=in2);
   by junction_id;
   if in2 then flag_junction_xscript_sim=1;
   else flag_junction_xscript_sim=0;
   if in1 then output;
run;

data sim_junctions_flagged2;
   set sim_junctions_flagged;
   if simulation="sim1" and test="test1" and flag_sim1_test1_xs=1 then flag_transcript_identified=1;
   else if simulation="sim2" and test="test1" and flag_sim2_test1_xs=1 then flag_transcript_identified=1;
   else if simulation="sim3" and test="test1" and flag_sim3_test1_xs=1 then flag_transcript_identified=1;
   else if simulation="sim1" and test="test2" and flag_sim1_test2_xs=1 then flag_transcript_identified=1;
   else if simulation="sim2" and test="test2" and flag_sim2_test2_xs=1 then flag_transcript_identified=1;
   else if simulation="sim3" and test="test2" and flag_sim3_test2_xs=1 then flag_transcript_identified=1;
   else flag_transcript_identified=0;
   keep simulation test junction_id flag_transcript_identified;
run;
 
/* Flag simulated junctions that are in the catalog */

data junc_in_cat;
   set event.event2star2pacbio_junc_table;
   where flag_in_catalog=1;
   keep junction_id;
run;

proc sort data=sim_junctions_flagged2;
   by junction_id;
proc sort data=junc_in_cat nodup;
   by junction_id;
run;

data sim_junctions_flagged3;
  merge sim_junctions_flagged2 (in=in1) junc_in_cat (in=in2);
  by junction_id;
  if in2 then flag_junction_in_catalog=1;
  else flag_junction_in_catalog=0;
  if in1 then output;
run;

/* Merge in APN counts for simulated junctions */

data junc_cat;
     set event.catalog_junctions_benchmark_sim;
     keep sample_id event_id apn;
     rename event_id=unique_junc_id;
run;

proc sort data=junc_cat;
   by unique_junc_id sample_id ;
proc transpose data=junc_cat out=junc_cat_sbys;
   by unique_junc_id;
   id sample_id;
   var apn;
run;

data seq2junc;
  set evspl.simulation_junc2uniq_seq;
  keep unique_junc_id junction_id;
run;

proc sort data=junc_cat_sbys;
  by unique_junc_id;
proc sort data=seq2junc nodup;
   by unique_junc_id;
run;

data junc_cat_sbys2;
  merge junc_cat_sbys (in=in1) seq2junc (in=in2);
  by unique_junc_id;
  if in1 and in2;
  drop unique_junc_id _NAME_;
run;

proc sort data=junc_cat_sbys2;
   by junction_id;
proc transpose data=junc_cat_sbys2 out=junc_cat2;
   by junction_id;
run;

data junc_cat3;
  length simulation $4;
  length test $5;
  set junc_cat2;
  simulation=scan(_NAME_,1,"_");
  test=scan(_NAME_,2,"_");
  drop _NAME_;
run;

proc sort data=junc_cat3;
   by simulation test junction_id;
proc sort data=sim_junctions_flagged3;
   by simulation test junction_id;
run;

data sim_junc_cat;
  merge sim_junctions_flagged3 (in=in1) junc_cat3 (in=in2);
  by simulation test junction_id;
  if in1 then flag_junction_sim_xscript=1; else flag_junction_sim_xscript=0;
  if in2 then flag_junc_in_catalog=1; else flag_junc_in_catalog=0;
  drop flag_junction_in_catalog;
run;

/* Merge in STAR counts for simulated junctions */

data junc_star;
  length simulation $4;
  length test $5;
  set event.star_junctions_benchmark_sim;
  if max_overhang lt 38 then delete;
  if strand=1 then strand_str="+"; else strand_str="-";
  donor_site=intron_start-1;
  acceptor_site=intron_stop;
  junction_id=catx(":",chr,donor_site,acceptor_site,strand_str);
  simulation=scan(sample_id,1,"_");
  test=scan(sample_id,2,"_");
  keep simulation test junction_id num_unique_mapped_reads;
  rename num_unique_mapped_reads=reads_mapped_star;
run;

proc sort data=junc_star;
   by simulation test junction_id;
proc sort data=sim_junc_cat;
   by simulation test junction_id;
run;

data sim_junc_cat2star;
  merge sim_junc_cat (in=in1) junc_star (in=in2) ;
  by simulation test junction_id;
  if in2 then flag_in_star=1; else flag_in_star=0;
  if not in1 then flag_junc_in_catalog=0;
run;

data sim_junc_cat2star2;
   set sim_junc_cat2star;
   if apn>0 and apn ne . then flag_events_detected=1;
   else flag_events_detected=0;
   if reads_mapped_star > 0 and reads_mapped_star ne . then flag_star_detected=1;
   else flag_star_detected=0;
   if flag_transcript_identified=. then flag_transcript_identified=0;
   if flag_junction_sim_xscript=. then flag_junction_sim_xscript=0;
run;

/* Count overlaps by sample */

proc sort data=sim_junc_cat2star2;
   by simulation test;
proc freq data=sim_junc_cat2star2 noprint;
   by simulation test;
   tables flag_transcript_identified*flag_junction_sim_xscript*flag_junc_in_catalog*flag_in_star*
          flag_events_detected*flag_star_detected / out=junc_counts_by_sample;
   tables flag_junction_sim_xscript*flag_events_detected*flag_star_detected / out=junc_counts_by_sample_simple;

run;

proc freq data=sim_junc_cat2star2 noprint;
   where flag_junc_in_catalog=1;
   by simulation test;
   tables flag_transcript_identified*flag_junction_sim_xscript*flag_junc_in_catalog*flag_in_star*
          flag_events_detected*flag_star_detected / out=junc_counts_by_sample2;
   tables flag_junction_sim_xscript*flag_events_detected*flag_star_detected / out=junc_counts_by_sample_simple2;
run;

proc print data=junc_counts_by_sample;
  where simulation="sim1" and test="test1";
run;

proc print data=junc_counts_by_sample_simple;
  where simulation="sim1" and test="test1";
run;

proc print data=junc_counts_by_sample2;
  where simulation="sim1" and test="test1";
run;

proc print data=junc_counts_by_sample_simple2;
  where simulation="sim1" and test="test1";
run;



/*
                                                           flag_
                       flag_junction_    flag_events_      star_
simulation    test       sim_xscript       detected      detected     COUNT

   sim1       test1           0                0             0       174677
   sim1       test1           0                0             1          170
   sim1       test1           0                1             0        20000
   sim1       test1           0                1             1           90
   sim1       test1           1                0             0       139171
   sim1       test1           1                0             1        10795
   sim1       test1           1                1             0        20195
   sim1       test1           1                1             1        86098


Or on just the junctions in the catalog:

                                    flag_
flag_junction_    flag_events_      star_
  sim_xscript       detected      detected     COUNT

       0                0             0       174677
       0                1             0        20000
       0                1             1           90
       1                0             0        28787
       1                0             1          429
       1                1             0        20195
       1                1             1        86098

*/


proc export data=junc_counts_by_sample
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_benchmarked_junctions_by_sample.csv"
     dbms=csv replace;
run;



proc export data=junc_counts_by_sample_simple
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_benchmarked_junctions_by_sample_simple.csv"
     dbms=csv replace;
run;




proc export data=junc_counts_by_sample2
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_benchmarked_junctions_by_sample_catalog_only.csv"
     dbms=csv replace;
run;



proc export data=junc_counts_by_sample_simple2
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_benchmarked_junctions_by_sample_simple_catalog_only.csv"
     dbms=csv replace;
run;


