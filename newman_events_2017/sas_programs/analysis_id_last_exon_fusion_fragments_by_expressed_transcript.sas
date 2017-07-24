ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* 3' variability check: check fragments of 3'-most exons matching to transcripts

For each transcript:
1. flag the fusion containing its last exon
2. flag the fragments of its last exon and number the order of fragments 5'->3' relative to orientation of transcript
3. flag last fragment
4. flag transcript if no detection in its last fusion
5. If last fusion detected, flag each fragment if not detected and cat together fragment ## detected and not detected

So, per transcript:

transcript
last fusion
last fragments
last fragment
last fusion dtct
# last fragments total
# last fragments dtct
cat fragments dtct (5>3)
cat fragments not dtct (5>3)

If all fragments are on, then no 3' variability
If some fragments off, then 3' variability

I also want to compare this to BLAST hits for PacBio later (need to think on how to do this)

I will probably later limit this to a subset of transcripts (e.g. ones with unique pieces)
but for now I will work on everything with at least one feature detected
*/

/* Get set of "expressed" transcripts */

data exp_xs;
   set event.feature_dtct_cnt_by_xscript_exp;
   where num_total_features_dtct > 0;
   keep transcript_id;
run;

data xs2exon;
   set mm10.mm10_exon2transcript;
run;

data exon_info;
   set mm10.mm10_exons_w_info;
   keep chrom start stop strand exon_id;
run;

proc sort data=xs2exon;
   by exon_id;
proc sort data=exon_info;
   by exon_id;
run;

data xs2exon_w_info;
  merge xs2exon (in=in1) exon_info (in=in2);
  by exon_id;
  if in1 and in2;
run;

proc sort data=exp_xs;
  by transcript_id;
proc sort data=xs2exon_w_info;
  by transcript_id;
run;

data exp_xs2exon_info;
  merge exp_xs (in=in1) xs2exon_w_info (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* strand check */

data xs2strand;
   set exp_xs2exon_info;
   keep transcript_id strand;
run;

proc sort data=xs2strand nodup;
   by transcript_id strand;
run; *101720, same number of xscripts as we start with;


proc sort data=exp_xs2exon_info;
   by transcript_id chrom start stop strand exon_id;
run;

data last_exon_per_xs;
   set exp_xs2exon_info;
   by transcript_id;
   if first.transcript_id and strand="-" then output;
   if last.transcript_id and strand="+" then output;
run;

/* Get fusion containing last exon
   I am only going to look at non-multigene fusions. If the 3' UTR of a transcript
   is in a multigene fusion, I will drop the transcript */

data fus2ex;
  set mm10.mm10_refseq_fusion_si_info_v2;
  if flag_multigene=0;
  keep fusion_id exon_id;
run;

proc sort data=fus2ex nodup;
  by exon_id fusion_id;
proc sort data=last_exon_per_xs;
  by exon_id;
run;

data last_fusion_per_xs last_exon_in_multigene;
  merge last_exon_per_xs (in=in1) fus2ex (in=in2);
  by exon_id;
  if in1 and in2 then output last_fusion_per_xs;
  else if in1 then output last_exon_in_multigene;
run;

*73535 of 73535 transcripts we can look at;

/* For each fusion, get the set of fragments. Flag fragments within fusion that are also assigned to transcript of interest */

data fus2xs;
   set last_fusion_per_xs;
   keep transcript_id fusion_id;
run;

proc sort data=fus2xs nodup;
   by fusion_id transcript_id;
proc freq data=fus2xs noprint;
   tables fusion_id / out=xs_per_fus;
proc sort data=xs_per_fus;
   by descending count;
proc print data=xs_per_fus(obs=1);
run; *60 max;

data cat_xscripts;
  array xs[60] $ 25;
  retain xs1-xs60;
  set fus2xs;
  by fusion_id;
  if first.fusion_id then do;
     call missing(of xs1-xs60);
     records = 0;
  end;
  records + 1;
  xs[records]=transcript_id;
  if last.fusion_id then output;
run;

data cat_xscripts2;
  set cat_xscripts;
  length transcript_cat $ 1000;
         transcript_cat= catx("|", OF xs1-xs60);
  keep fusion_id transcript_cat;
  run;

data fus2frag;
   set mm10.mm10_exon_fragment_flagged;
   keep fragment_id fusion_id transcript_id exon_id fragment_Start fragment_end;
   rename exon_id=exons_of_frag transcript_id=transcripts_of_frag;
run;

proc sort data=fus2frag nodup;
   by  fusion_id fragment_id;
proc sort data=cat_xscripts2;
   by fusion_id;
run;

data last_frags_per_xs;
  merge cat_xscripts2 (in=in1) fus2frag (in=in2);
  by fusion_id;
  if in1 and in2;
run;

/* uncat transcripts */

data last_frags_per_xs2;
   length transcript_id $25.;
   set last_frags_per_xs;
   do i=1 by 1 while(scan(transcript_cat,i,"|") ^= "");
      transcript_id=scan(transcript_cat,i,"|");
      output;
      end;
   drop transcript_cat i;
run;

proc sort data=last_frags_per_xs2;
   by transcript_id fusion_id;
proc sort data=last_fusion_per_xs;
   by transcript_id fusion_id;
run;

data last_frags_per_xs3;
   merge last_fusion_per_xs (in=in1) last_frags_per_xs2 (in=in2);
   by transcript_id fusion_id;
   if in1 and in2;
run;

/* Flag if fragment is assigned to transcript and last exon of transcript */

data flag_frag;
   set last_frags_per_xs3;
   if index(exons_of_frag,strip(exon_id)) > 0 then flag_frag_of_last_exon=1;
   else flag_frag_of_last_exon=0;
   if index(transcripts_of_frag,strip(transcript_id)) > 0 then flag_frag_of_xs=1;
   else flag_frag_of_xs=0;
run;

proc freq data=flag_frag;
   tables flag_frag_of_xs*flag_frag_of_last_exon;
run;

/*
   flag_frag_of_xs
             flag_frag_of_last_exon

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  28235 |     46 |  28281
            |  19.24 |   0.03 |  19.27
            |  99.84 |   0.16 |
            |  84.47 |   0.04 |
   ---------+--------+--------+
          1 |   5193 | 113288 | 118481
            |   3.54 |  77.19 |  80.73
            |   4.38 |  95.62 |
            |  15.53 |  99.96 |
   ---------+--------+--------+
   Total       33428   113334   146762
               22.78    77.22   100.00


*/

/* Check that all transcripts have at least one fragment */

proc sort data=flag_frag;
   by transcript_id;
proc means data=flag_frag noprint;
   by transcript_id;
   var flag_frag_of_last_exon flag_frag_of_xs;
   output out=xs_frag_check max=;
run;

proc freq data=xs_frag_check;
   tables flag_frag_of_xs*flag_frag_of_last_exon;
run;

/*
 flag_frag_of_xs
           flag_frag_of_last_exon

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        1 |  73535 |  73535
          | 100.00 | 100.00
          | 100.00 |
          | 100.00 |
 ---------+--------+
 Total       73535    73535
            100.00   100.00


*/

/* Make permenant */

data event.last_frag_exon_fus_by_xscript;
    set flag_frag;
run;



