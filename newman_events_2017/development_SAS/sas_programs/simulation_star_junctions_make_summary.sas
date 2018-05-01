
/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';


/* For each test simulation, I want to count the number of annotated junctions, unannotated junction, novel junctions
   not detected, detected in 1 sample, 2 samples, all 3 samples, as well as their proportions */

data star_junc_w_annot;
  set  event.simul_star_junc_w_annot;
  length junction_category $15.;
  if flag_junc_annotated_cat=1 then junction_category="annotated";
  else if flag_junc_annotated_cat=0 then junction_category="unannotated";
  else junction_category="novel";
run;

%macro cntJunc(test);
/* Count total number of junctions by type */
proc freq data=star_junc_w_annot noprint;
  tables junction_category / out=all_junc_cnt;
run;

/* Count total number of junctions without coverage in any of the test1 samples */
proc freq data=star_junc_w_annot noprint;
  where sim1_&test.=0 and sim2_&test.=0 and sim3_&test.=0;
  tables junction_category / out=null_junc_cnt;
run;

/* Count total number of junctions with coverage in the test1 samples, but not considered "detected" in any */
proc freq data=star_junc_w_annot noprint;
  where sum(sim1_&test.,sim2_&test.,sim3_&test.) > 0 and flag_depth_sim1_&test._ge5=0 
        and flag_depth_sim2_&test._ge5=0 and flag_depth_sim3_&test._ge5=0 ;
  tables junction_category / out=undtct_junc_cnt;
run;

/* Count total number of junctions detected in only one sample  */
proc freq data=star_junc_w_annot noprint;
  where sum(flag_depth_sim1_&test._ge5, flag_depth_sim2_&test._ge5, flag_depth_sim3_&test._ge5) = 1 ;
  tables junction_category / out=n1_junc_cnt;
run;


/* Count total number of junctions detected in two samples  */
proc freq data=star_junc_w_annot noprint;
  where sum(flag_depth_sim1_&test._ge5, flag_depth_sim2_&test._ge5, flag_depth_sim3_&test._ge5) = 2 ;
  tables junction_category / out=n2_junc_cnt;
run;

/* Count total number of junctions detected in all three samples  */
proc freq data=star_junc_w_annot noprint;
  where sum(flag_depth_sim1_&test._ge5, flag_depth_sim2_&test._ge5, flag_depth_sim3_&test._ge5) = 3 ;
  tables junction_category / out=n3_junc_cnt;
run;

data all_junc_cnt2;
  set all_junc_cnt;
  length detection $60.;
  detection="Total junctions in test set";
run;

data null_junc_cnt2;
  set null_junc_cnt;
  length detection $60.;
  detection="Junctions not in test set";
run;

data undtct_junc_cnt2;
  set undtct_junc_cnt;
  length detection $60.;
  detection="Junctions undetected";
run;

data n1_junc_cnt2;
  set n1_junc_cnt;
  length detection $60.;
  detection="Junctions detected in 1 sample";
run;

data n2_junc_cnt2;
  set n2_junc_cnt;
  length detection $60.;
  detection="Junctions detected in 2 samples";
run;

data n3_junc_cnt2;
  set n3_junc_cnt;
  length detection $60.;
  detection="Junctions detected in 3 samples";
run;

proc transpose data=all_junc_cnt2 out=all_junc_cnt_sbys(drop=_NAME_ _LABEL_);
   by detection;
   var count;
   id junction_category;
run;

proc transpose data=null_junc_cnt2 out=null_junc_cnt_sbys(drop=_NAME_ _LABEL_);
   by detection;
   var count;
   id junction_category;
run;

proc transpose data=undtct_junc_cnt2 out=undtct_junc_cnt_sbys(drop=_NAME_ _LABEL_);
   by detection;
   var count;
   id junction_category;
run;

proc transpose data=n1_junc_cnt2 out=n1_junc_cnt_sbys(drop=_NAME_ _LABEL_);
   by detection;
   var count;
   id junction_category;
run;
proc transpose data=n2_junc_cnt2 out=n2_junc_cnt_sbys(drop=_NAME_ _LABEL_);
   by detection;
   var count;
   id junction_category;
run;

proc transpose data=n3_junc_cnt2 out=n3_junc_cnt_sbys(drop=_NAME_ _LABEL_);
   by detection;
   var count;
   id junction_category;
run;

data &test._junctions;
   set all_junc_cnt_sbys null_junc_cnt_sbys undtct_junc_cnt_sbys 
   n1_junc_cnt_sbys n2_junc_cnt_sbys n3_junc_cnt_sbys;
   total=annotated + novel + unannotated;
   perc_annotated=annotated/total * 100;
   perc_novel=novel/total * 100;
   perc_unannotated=unannotated/total * 100;
run;

proc print data=&test._junctions;
run;

%mend;

%cntJunc(test1);
%cntJunc(test2);

/*  TEST 1:

Obs    detection                                 annotated     novel

 1     Total junctions in set                      181249     109080
 2     Junctions not in Test 1 set                  20661      54315
 3     Test 1 junctions undetected                  38300      28822
 4     Test 1 junctions detected in 1 sample        43601      23090
 5     Test 1 junctions detected in 2 samples       44823       2628
 6     Test 1 junctions detected in 3 samples       33864        225

                                  perc_       perc_        perc_
Obs    unannotated     total    annotated     novel     unannotated

 1        195501      485830     37.3071     22.4523      40.2406
 2         83264      158240     13.0567     34.3244      52.6188
 3         73257      140379     27.2833     20.5316      52.1852
 4         33162       99853     43.6652     23.1240      33.2108
 5          4934       52385     85.5646      5.0167       9.4187
 6           884       34973     96.8290      0.6434       2.5277

TEST 2:

 Obs    detection                          annotated     novel

  1     Total junctions in test set          181249     109080
  2     Junctions not in test set             23104      52993
  3     Junctions undetected                  40533      30831
  4     Junctions detected in 1 sample        42799      22829
  5     Junctions detected in 2 samples       43594       2238
  6     Junctions detected in 3 samples       31219        189

                                   perc_       perc_        perc_
 Obs    unannotated     total    annotated     novel     unannotated

  1        195501      485830     37.3071     22.4523      40.2406
  2         71219      147316     15.6833     35.9723      48.3444
  3         70391      141755     28.5937     21.7495      49.6568
  4         46031      111659     38.3301     20.4453      41.2246
  5          6926       52758     82.6301      4.2420      13.1279
  6           934       32342     96.5277      0.5844       2.8879

Assume "high confidence" junctions are those detected at depth >= 5 unique reads in at least 2 of 3 samples,
then most detected junctions are annotated junctions, then unannotated junctions, then novel


Should test this on the mouse neural samples...
*/


