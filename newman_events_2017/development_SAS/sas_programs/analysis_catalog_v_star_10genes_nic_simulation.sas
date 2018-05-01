/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the 10-gene simulation, import the STAR junction results and EA junction coverage counts and check
   to see which method captures the junctions of these genes */

/* Junctions from genes simulated */

data genes_sim;
  set event.genes_w_nic_junction_10genes;
run;


data junc2check;
   set evspl.splicing_events_annot_refseq;
   length junction_id $50.;
   where flag_intron_retention=0;
   junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
   keep gene_id event_id junction_id flag_junction_annotated;
run;

proc sort data=genes_sim;
  by gene_id;
proc sort data=junc2check;
  by gene_id;
run;

data junc2check2;
  merge genes_sim (in=in1) junc2check (in=in2);
  by gene_id;
  if in1 and in2 ;
run;

data event2uniq;
  set eventloc.unique_junction2event_mm10;
  keep event_id junction_id seq_name;
run;

data cat_junc;
   set event.catalog_junctions_10genes_nic;
   keep sample_id event_id apn region_depth reads_in_region;
   rename event_id=seq_name apn=event_apn region_depth=event_depth reads_in_region=event_reads_in_region;
run;

proc sort data=cat_junc;
  by seq_name;
proc sort data=event2uniq;
  by seq_name;
run;

data cat_junc2;
  merge event2uniq (in=in1) cat_junc (in=in2);
  by seq_name;
  if in1 and in2;
run;

proc sort data=cat_junc2;
   by junction_id event_id;
proc sort data=junc2check2;
   by junction_id event_id;
run;

data cat_junc3;
   merge junc2check2 (in=in1) cat_junc2 (in=in2);
   by junction_id event_id;
   if in2;
  if gene_id ^= "" then flag_gene_sim=1; else flag_gene_sim=0;
run;


data star_junc;
  set event.star_junctions_10genes_nic;
  length junction_id $50.;
  if max_overhang lt 16 then delete;
  if strand=1 then strand_str="+"; else strand_str="-";
  donor_site=intron_start-1;
  acceptor_site=intron_stop;
  junction_id=catx(":",chr,donor_site,acceptor_site,strand_str);
  keep sample_id junction_id num_unique_mapped_reads;
  rename num_unique_mapped_reads=reads_mapped_star;
run;

proc sort data=star_junc;
   by sample_id junction_id;
proc sort data=cat_junc3;
   by sample_id junction_id;
run;

data star2cat_10genes;
  merge cat_junc3 (in=in1) star_junc (in=in2);
  by sample_id junction_id;
  if in1 then flag_in_catalog=1; else flag_in_catalog=0;
  if in2 then flag_in_star=1; else flag_in_star=0;
run;

data flag_on;
  set star2cat_10genes;
  if event_apn=0 or event_apn=. then flag_event_on=0;
  else flaG_event_on=1;
  if reads_mapped_star=0 or reads_mapped_star=. then flag_star_on=0;
  else flaG_star_on=1;
run;

data juncs_sim;
   set junc2check2;
   if flag_junction_annotated=1 then output;
   else if junction_id in ('chr6:147719923:147727642:-','chr9:48456829:48463053:+','chr6:119921643:119922471:+',
                           'chr5:53466228:53600832:+','chr16:36875263:36885010:+','chr3:88736088:88748781:-',
                           'chr5:129109661:129128512:+', 'chr7:4631022:4631122:-','chr7:132837207:132850590:-',
                           'chr13:106823051:106836211:-') then output;
   keep junction_id;
run;

proc sort data=flag_on;
  by junction_id;
proc sort data=juncs_sim;
  by junction_id;
run;

data flag_on2;
   merge flag_on (in=in1) juncs_sim (in=in2);
   by junction_id;
   if in2 then flag_junction_simulated=1;
   else flag_junction_simulated=0;
   if in1 then output;
run;


proc sort data=flag_on2;
  by sample_id;
proc freq data=flag_on2 noprint;
  by sample_id;
  tables flag_in_catalog*flag_in_star*flag_event_on*flag_star_on*flag_junction_simulated*flag_junction_annotated
         / out=junc_counts_by_sample;
run;

proc print data=junc_counts_by_sample;
where sample_id="sample_01";
run;

/* For sample 1:
           flag_in_ flag_in_   flag_   flag_  flag_junction_ flag_junction_
 sample_id  catalog   star   event_on star_on    simulated      annotated   COUNT PERCENT

 sample_01     1        0        0       0           0              .         17    .
 sample_01     1        0        1       0           0              .         28    .

 sample_01     1        0        1       0           0              0          1   0.5882
 sample_01     1        0        1       0           1              0          4   2.3529
 sample_01     1        0        1       0           1              1          2   1.1765

 sample_01     1        1        1       1           0              0          9   5.2941
 sample_01     1        1        1       1           1              0          5   2.9412
 sample_01     1        1        1       1           1              1        149  87.6471

In terms of detection:
Events: gets all junctions from simulated transcripts, plus 9 NIC junctions.
        It also finds 10 additional (NIC) junction
        Also finds evidence of an additional 45 border junctions
STAR: misses 2 annotated junctions, and 4 of the NIC junctions
        It also finds 1 additional (NIC) junction

*/

proc export data=junc_counts_by_sample
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_10genes_nic_junctions_by_sample.csv"
     dbms=csv replace;
run;





