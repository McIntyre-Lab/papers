/* Libraries */

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
ods listing; ods html close;

/* Stack together the iReckon results for each sample.
   I want to extract the exon sequences only, so I can
   drop the transcript entries */

/* first I want to count how many exons are in common between NSC1 and NSC2 */

data nsc1_exons;
   set event.iReckon_nsc1_output;
   where feature_type="exon";
   keep chr start stop strand;
run;

data nsc2_exons;
   set event.iReckon_nsc2_output;
   where feature_type="exon";
   keep chr start stop strand;
run;

proc sort data=nsc1_exons nodup;
   by chr start stop strand;
proc sort data=nsc2_exons nodup;
   by chr start stop strand;
run;
 
data nsc1v2_exons;
  merge nsc1_exons (in=in1) nsc2_exons (in=in2);
  by chr start stop strand;
  if in1 then flag_nsc1=1; else flag_nsc1=0;
  if in2 then flag_nsc2=1; else flag_nsc2=0;
run;

proc freq data=nsc1v2_exons;
  tables flag_nsc1*flag_nsc2;
run;

/*
    Table of flag_nsc1 by flag_nsc2

  flag_nsc1     flag_nsc2

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   9123 |   9123
           |   0.00 |  33.47 |  33.47
           |   0.00 | 100.00 |
           |   0.00 |  60.78 |
  ---------+--------+--------+
         1 |  12246 |   5888 |  18134
           |  44.93 |  21.60 |  66.53
           |  67.53 |  32.47 |
           | 100.00 |  39.22 |
  ---------+--------+--------+
  Total       12246    15011    27257
              44.93    55.07   100.00

So most distinct exons are unique to each sample
*/

data nsc_exons;
  set event.iReckon_nsc1_output (in=in1) event.iReckon_nsc2_output (in=in2) ;
  length sample_id $4.;
  where feature_type="exon";
  if in1 then sample_id="NSC1";
  else if in2 then sample_id="NSC2";
  keep chr start stop strand gene_id transcript_id exon_number sample_id;
run;

data xs_gene_new_ids;
   set nsc_exons;
   length gene_id2 $20.;
   length transcript_id2 $60.;
   length exon_id $60.;
   length exon_coord_id $30.;
   gene_id2=catx(":",sample_id,gene_id);
   transcript_id2=catx(":",sample_id,gene_id,transcript_id);
   exon_id=catx(":",sample_id,gene_id,transcript_id,exon_number);
   exon_coord_id=catx(":",chr,start,stop,strand);
   keep gene_id2 transcript_id2 exon_id exon_coord_id chr start stop strand;
   rename gene_id2=gene_id;
   rename transcript_id2=transcript_id;
run;

proc sort data=xs_gene_new_ids nodup ;
   by gene_id transcript_id exon_id exon_coord_id;
run;

/* This is my main "exon index" that I will use later. Make this permenant for now.
   Then I need to export the list of distinct exons for extracting sequences */

data event.ireckon_nsc_exon2xscript2gene;
  set xs_gene_new_ids;
  drop chr start stop strand;
run;

data distinct_exons;
  set xs_gene_new_ids;
  keep exon_coord_id  chr start stop strand;
run;

data distinct_exons_bed;
  retain chr start stop  exon_coord_id score strand;
  set  distinct_exons ;
  score=".";
  if strand="" then strand="+"; * if strand is blank, assume sense strand;
run; 

/* check that coordinates are around the right way */

proc sort data=distinct_exons_bed nodup;
  by chr start stop exon_coord_id score strand;
run; *27257 exons total to check;


data check_pos;
  set distinct_exons_bed;
  if start < stop then flag_pos_ok=1;
  else flag_pos_ok=0;
run;

proc freq data=check_pos;
  tables flag_pos_ok;
run;

proc print data=check_pos;
  where flaG_pos_ok=0;
run;

/*

                                              Cumulative    Cumulative
      flag_pos_ok    Frequency     Percent     Frequency      Percent
      ----------------------------------------------------------------
                0          10        0.04            10         0.04
                1       27247       99.96         27257       100.00

Only 10 of these. And they are all 1 bp long, so I will drop them for now
since we can't really say where they will align to

I might check the coordinates of these with RefSeq exons to see if they are annotated
*/

data distinct_exons_bed2;
  set check_pos;
  where flag_pos_ok=1;
  start=start-1;
  drop flag_pos_ok;
run;


proc export data=distinct_exons_bed2
     outfile="!MCLAB/event_analysis/analysis_output/ireckon_nsc_distinct_exons.bed"
     dbms=tab replace;
     putnames=no;
run;

