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

I also want to compare this to BLAST hits for PacBio later (need to think on how to do this)

I will probably later limit this to a subset of transcripts (e.g. ones with unique pieces)
but for now I will work on everything with at least one feature detected
*/


/* Fragment counts:
1. total fragments in region
2. total detected fragments in region
3. total fragments in region that are also in xscript
4. total detected fragments in region that are also in xscript
5. Cat together fragment orders:
	a. fragments in xscript 
	b. detected fragments in xscript
	c. fragments not in xscript
	d. detected fragments not in xscript

I also want to check that the set of detected fragments are sequential
*/

data frag2xs_frag_fus_on;
   set event.xscripts_last_fus_frag_onflags;
   where flag_fusion_on=1;
run;

/* Count total fragments in region */

proc sort data=frag2xs_frag_fus_on;
   by transcript_id ;
proc freq data=frag2xs_frag_fus_on noprint;
   tables transcript_id / out = frags_per_region;
run;

/* Count detected fragments in region */

proc freq data=frag2xs_frag_fus_on noprint;
   where flag_fragment_on=1;
   tables transcript_id / out = dtct_frags_per_region;
run;

/* Count total fragments in region that are also assigned to transcript */

proc freq data=frag2xs_frag_fus_on noprint;
   where flag_frag_of_xs=1;
   tables transcript_id / out = frags_per_xs;
run;

/* Count detected fragments in region that are also assigned to transcript */

proc freq data=frag2xs_frag_fus_on noprint;
   where flag_frag_of_xs=1 and flag_fragment_on=1;
   tables transcript_id / out = dtct_frags_per_xs;
run;

/* Cat together fragment numbers */

* fragments in xscript;
data xs2frag;
   set frag2xs_frag_fus_on;
   where flag_frag_of_xs=1;
   keep transcript_id fragment_id frag_num;
run;

proc sort data=xs2frag nodup;
   by transcript_id frag_num;
proc freq data=xs2frag noprint;
   tables transcript_id / out=frag_per_xs;
proc sort data=frag_per_xs;
   by descending count;
proc print data=frag_per_xs(obs=1);
run; *41 max;

data cat_frag_num;
  array frag[41] $ 2;
  retain frag1-frag41;
  set xs2frag;
  by transcript_id;
  if first.transcript_id then do;
     call missing(of frag1-frag41);
     records = 0;
  end;
  records + 1;
  frag[records]=frag_num;
  if last.transcript_id then output;
run;

data cat_frags_per_xs;
  set cat_frag_num;
  length xscript_fragments_cat $ 150;
         xscript_fragments_cat= catx("|", OF frag1-frag41);
  keep transcript_id xscript_fragments_cat;
  run;


* detected fragments in xscript;
data xs2frag;
   set frag2xs_frag_fus_on;
   where flag_frag_of_xs=1 and flag_fragment_on=1;
   keep transcript_id fragment_id frag_num;
run;

proc sort data=xs2frag nodup;
   by transcript_id frag_num;
proc freq data=xs2frag noprint;
   tables transcript_id / out=frag_per_xs;
proc sort data=frag_per_xs;
   by descending count;
proc print data=frag_per_xs(obs=1);
run; *18 max;

data cat_frag_num;
  array frag[18] $ 2;
  retain frag1-frag18;
  set xs2frag;
  by transcript_id;
  if first.transcript_id then do;
     call missing(of frag1-frag18);
     records = 0;
  end;
  records + 1;
  frag[records]=frag_num;
  if last.transcript_id then output;
run;

data cat_dtct_frags_per_xs;
  set cat_frag_num;
  length dtct_xscript_fragments_cat $ 150;
         dtct_xscript_fragments_cat= catx("|", OF frag1-frag18);
  keep transcript_id dtct_xscript_fragments_cat;
  run;


* fragments not in xscript;
data xs2frag;
   set frag2xs_frag_fus_on;
   where flag_frag_of_xs=0;
   keep transcript_id fragment_id frag_num;
run;


proc sort data=xs2frag nodup;
   by transcript_id frag_num;
proc freq data=xs2frag noprint;
   tables transcript_id / out=frag_per_xs;
proc sort data=frag_per_xs;
   by descending count;
proc print data=frag_per_xs(obs=1);
run; *25 max;

data cat_frag_num;
  array frag[25] $ 2;
  retain frag1-frag25;
  set xs2frag;
  by transcript_id;
  if first.transcript_id then do;
     call missing(of frag1-frag25);
     records = 0;
  end;
  records + 1;
  frag[records]=frag_num;
  if last.transcript_id then output;
run;

data cat_frags_not_in_xs;
  set cat_frag_num;
  length non_xscript_fragments_cat $ 150;
         non_xscript_fragments_cat= catx("|", OF frag1-frag25);
  keep transcript_id non_xscript_fragments_cat;
  run;

* detected fragments not in xscript;
data xs2frag;
   set frag2xs_frag_fus_on;
   where flag_frag_of_xs=0 and flag_fragment_on=1;
   keep transcript_id fragment_id frag_num;
run;


proc sort data=xs2frag nodup;
   by transcript_id frag_num;
proc freq data=xs2frag noprint;
   tables transcript_id / out=frag_per_xs;
proc sort data=frag_per_xs;
   by descending count;
proc print data=frag_per_xs(obs=1);
run; *14 max;

data cat_frag_num;
  array frag[14] $ 2;
  retain frag1-frag14;
  set xs2frag;
  by transcript_id;
  if first.transcript_id then do;
     call missing(of frag1-frag14);
     records = 0;
  end;
  records + 1;
  frag[records]=frag_num;
  if last.transcript_id then output;
run;

data cat_dtct_frags_not_in_xs;
  set cat_frag_num;
  length dtct_non_xscript_fragments_cat $ 150;
         dtct_non_xscript_fragments_cat= catx("|", OF frag1-frag14);
  keep transcript_id dtct_non_xscript_fragments_cat;
  run;

/* Merge all together and make permenant */

proc sort data=frags_per_region (rename=(count=frags_per_region));
  by transcript_id;
proc sort data=dtct_frags_per_region (rename=(count=dtct_frags_per_region));
  by transcript_id;
proc sort data=frags_per_xs (rename=(count=frags_per_xscript));
  by transcript_id;
proc sort data=dtct_frags_per_xs (rename=(count=dtct_frags_per_xscript));
  by transcript_id;
proc sort data=cat_frags_per_xs;
  by transcript_id;
proc sort data=cat_dtct_frags_per_xs;
  by transcript_id;
proc sort data=cat_frags_not_in_xs;
  by transcript_id;
proc sort data=cat_dtct_frags_not_in_xs;
  by transcript_id;
run;

data xscript_w_frag_lists;
  merge frags_per_region dtct_frags_per_region frags_per_xs dtct_frags_per_xs
        cat_frags_per_xs cat_dtct_frags_per_xs  cat_frags_not_in_xs cat_dtct_frags_not_in_xs;
  by transcript_id;
run;

data xscript_w_frag_lists2;
   set xscript_w_frag_lists;
   array change _numeric_;
            do over change;
            if change=. then change=0;
            end;
   drop PERCENT;
run;

data flag_3prime_var;
   set xscript_w_frag_lists2;
   if dtct_frags_per_xscript = frags_per_xscript then flag_possible_3prime_var=0;
   else flag_possible_3prime_var=1;
run;

data xscripts_last_frag_flag;
   set frag2xs_frag_fus_on;
   where flag_frag_of_xs=1 and flag_frag_last_in_xscript=1;
   keep transcript_id flag_fragment_on;
   rename flag_fragment_on=flag_last_frag_in_xscript_on;
run;

proc sort data=flag_3prime_var;
   by transcript_id;
proc sort data=xscripts_last_frag_flag;
   by transcript_id;
run;

data event.xscripts_3prime_fragments;
   merge flag_3prime_var (in=in1) xscripts_last_frag_flag (in=in2);
   by transcript_id;
   if in1 and in2;
run;

proc freq data=event.xscripts_3prime_fragments;
   tables flag_possible_3prime_var
          flag_last_frag_in_xscript_on
          flag_possible_3prime_var*flag_last_frag_in_xscript_on  ;
run;


/*

     flag_possible_                             Cumulative    Cumulative
         3prime_var    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0       57718       90.46         57718        90.46
                  1        6089        9.54         63807       100.00


    flag_last_frag_                             Cumulative    Cumulative
      in_xscript_on    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        3064        4.80          3064         4.80
                  1       60743       95.20         63807       100.00


 flag_possible_3prime_var
           flag_last_frag_in_xscript_on

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |  57718 |  57718
          |   0.00 |  90.46 |  90.46
          |   0.00 | 100.00 |
          |   0.00 |  95.02 |
 ---------+--------+--------+
        1 |   3064 |   3025 |   6089
          |   4.80 |   4.74 |   9.54
          |  50.32 |  49.68 |
          | 100.00 |   4.98 |
 ---------+--------+--------+
 Total        3064    60743    63807
              4.80    95.20   100.00


3064 transcripts with 3' variability, and the last fragment of the transcript is not detected
*/

