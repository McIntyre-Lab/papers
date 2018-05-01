/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Of the analyzable PacBio transcripts, what is the detection rate in terms of novel junctions/isoforms
   for STAR, Event Analysis, and iReckon?

   Junction-based:

   (1) Identify novel PB junctions, and their associated transcripts
   (2) For STAR/EA, how many novel junctions are detected and supported?
   (3) Collapse (2) back to (a) transcripts and (b) genes:
		(a) How many novel transcripts are detected by STAR and EA?
		(b) How many genes with novel transcripts are detected by STAR and EA?


   Transcript-based:
   (1) Identify the set of novel PacBio transcripts
   (2) What is the overlap between iReckon transcripts and novel PacBio?
		-> NSC1, NSC2, NSC1 ∪ NSC2, NSC1 ∩ NSC2
		-> Cross-tabulation
	(note this will be more-or-less the same for the general iReckon vs PacBio comparison
     so code here can be used again)
*/

/*** JUNCTION-BASED ***/

/* (1) Identify novel PB junctions
   I am not including annotated junctions that are present in the novel transcripts, as we cannot
   disentangle these from annotated transcripts, even IF the novel transcript is the only
   transcript present.

   I might need to flag these transcripts: ie. transcripts without NNC or NIC junctions . Should
   be zero, unless it retains an intron */

data pb_xs;
   set event.pacbio_transcripts;
run;

proc freq data=pb_xs;
  tables pb_type;
run;

/* Overall PB transcripts by type (including genes with multigene pieces):

pb_type                    Frequency     Percent
-------------------------------------------------
antisense                       177        1.10
full-splice_match              7938       49.29
fusion                           54        0.34
genic                            93        0.58
genic_intron                    299        1.86
incomplete-splice_match        1742       10.82
intergenic                       61        0.38
novel_in_catalog               3045       18.91
novel_not_in_catalog           2695       16.73

*/

data pb_xs_nomult;
   set event.pacbio_transcripts_nomulti;
run;

proc freq data=pb_xs_nomult;
  tables pb_type;
run;

/* Overall PB transcripts by type (including genes with multigene pieces):

pb_type                    Frequency     Percent
---------------------------------------------------
antisense                        10        0.10
full-splice_match              5728       54.56
fusion                           32        0.30
genic                            14        0.13
genic_intron                     54        0.51
incomplete-splice_match        1213       11.55
intergenic                        1        0.01
novel_in_catalog               1730       16.48
novel_not_in_catalog           1716       16.35

*/

data pb_novel;
   set event.pacbio_transcripts;
   where pb_type ^? "match";
   keep pacbio_id;
run;

data pb_novel_nomulti;
   set event.pacbio_transcripts_nomulti;
   where pb_type ^? "match";
   keep pacbio_id;
run;

proc sort data=pb_novel;
   by pacbio_id;
proc sort data=pb_novel_nomulti;
   by pacbio_id;
run;

data pb_novel_flag_multi;
  merge pb_novel (in=in1) pb_novel_nomulti (in=in2);
  by pacbio_id;
  if in2 then flag_multigene=0; else flag_multigene=1;
run;

/* Flag which of these have junctions */

data pb_junc;
  set evspl.splicing_events_annotations;
  length junction_id $50.;
  length pacbio_id $15.;
  junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      pacbio_id=scan(transcript_id,i,"|");
      output;
      end;
  keep junction_id pacbio_id;
run;

proc sort data=pb_junc nodup;
   by junction_id pacbio_id;
run;

data pb_w_junc;
  set pb_junc;
  keep pacbio_id;
run;

proc sort data=pb_w_junc nodup;
   by pacbio_id;  *14657 PB transcripts with junctions;
proc sort data=pb_novel_flag_multi;
   by pacbio_id;
run;

data pb_novel_flag_junc;
  merge pb_novel_flag_multi (in=in1) pb_w_junc (in=in2);
  by pacbio_id;
  if in2 then flag_has_junction=1; else flag_has_junction=0;
  if in1 then output;
run;

proc freq data=pb_novel_flag_junc;
  tables flag_multigene*flag_has_junction;
run;

/*
 flag_multigene
           flag_has_junction

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    172 |   3385 |   3557
          |   2.68 |  52.69 |  55.37
          |   4.84 |  95.16 |
          |  22.16 |  59.93 |
 ---------+--------+--------+
        1 |    604 |   2263 |   2867
          |   9.40 |  35.23 |  44.63
          |  21.07 |  78.93 |
          |  77.84 |  40.07 |
 ---------+--------+--------+
 Total         776     5648     6424
             12.08    87.92   100.00



3385 of 3557 novel PB transcripts not in multigene regions have junctions
5648 of 6424 novel PB transcripts (total) have junctions

172-776 transcripts will be missed at the junction regardless of method!
*/

/* Flag novel transcripts that have novel junctions (NIC or NNC) */

data nnc_junc nic_junc;
   set event.event2star2pacbio_junc_table;
   where flag_in_pacbio=1 ;
   if junction_id="" then junction_id=catx(":",chr,donor_stop,acceptor_start,strand);
   if index(junction_type,"NNC") > 0 then output nnc_junc;
   else if index(junction_type,"NIC") > 0 then output nic_junc;
   keep junction_id junction_type;
run;

proc sort data=nnc_junc nodup;
   by junction_id;
proc sort data=nic_junc nodup;
   by junction_id;
proc sort data=pb_junc;
   by junction_id;
run;

data pb_nnc_junc;
   merge pb_junc (in=in1) nnc_junc (in=in2);
   by junction_id;
   if in1 and in2;
run;

data pb_nic_junc;
   merge pb_junc (in=in1) nic_junc (in=in2);
   by junction_id;
   if in1 and in2;
run;


data pb_novel_w_nnc_junc;
   set pb_nnc_junc;
   keep pacbio_id;
run;

data pb_novel_w_nic_junc;
   set pb_nic_junc;
   keep pacbio_id;
run;

proc sort data=pb_novel_w_nnc_junc nodup;
   by pacbio_id; *3213 transcripts;
proc sort data=pb_novel_w_nic_junc nodup;
   by pacbio_id; *579 transcripts;
proc sort data=pb_novel_flag_junc;
   by pacbio_id;
run;

data pb_novel_flag_novel_junc;
  merge pb_novel_flag_junc (in=in1) pb_novel_w_nnc_junc (in=in2) pb_novel_w_nic_junc (in=in3);
  by pacbio_id;
  if in2 or in3 then flag_has_novel_junction=1; else flag_has_novel_junction=0;
  if in2 then flag_has_nnc_junction=1; else flag_has_nnc_junction=0;
  if in3 then flag_has_nic_junction=1; else flag_has_nic_junction=0;
  if in1;
run;

proc freq data=pb_novel_flag_novel_junc noprint;
  tables flag_multigene*flag_has_junction*flag_has_novel_junction / out=pb_junc_count;
run;

proc print data=pb_junc_count;
run;

/*
   flag_      flag_has_      novel_
 multigene     junction     junction    COUNT

     0            0            0          172
     0            1            0         1149
     0            1            1         2236
     1            0            0          604
     1            1            0          915
     1            1            1         1348

172+604 novel isoforms without junctions
1149+915 novel isoforms made from completely known junctions
2236 + 1348 novel isoforms with at least one novel junction

By type?
*/

proc freq data=pb_novel_flag_novel_junc noprint;
  tables flag_multigene*flag_has_junction*flag_has_nnc_junction*flag_has_nic_junction / out=pb_junc_count;
run;

proc print data=pb_junc_count;
run;

/*
    flag_      flag_has_      flag_has_       flag_has_
  multigene     junction    nnc_junction    nic_junction    COUNT
NO JUNCTIONS:
      0            0              0               0           172
      1            0              0               0           604

NO NOVEL JUNCTIONS:
      0            1              0               0          1149
      1            1              0               0           915

NNC JUNCTIONS:
      0            1              1               0          1884
      1            1              1               0          1154

NIC JUNCTIONS
      0            1              0               1           326
      1            1              0               1           161

NNC+NIC JUNCTIONS:
      0            1              1               1            26
      1            1              1               1            33

172 + 604 novel isoforms with no junctions
1149 + 915 novel isoforms with no novel junctions
(172+604+1149+915 isoforms we can't identify from junctions)
1884+1154 novel isoforms with NNC junctions "only" (may have annotated junctions)
326+161 novel isoforms with NIC junctions "only" (may have annotated junctions)
26+33 novel isoforms with NIC and NNC junctions
*/

/* Get novel junctions that are detected/have support in EA/STAR */

data junc_flags;
   set event.event2star2pacbio_junc_table;
   where flag_in_pacbio=1 and (junction_type ? "NNC" or junction_type ? "NIC");
   if junction_id="" then junction_id=catx(":",chr,donor_stop,acceptor_start,strand);
   keep junction_id junction_type flag_events_detected flag_events_NSC_apn_gt0
        flag_events_NSC_apn_ge2 flag_events_NSC_apn_ge5 flag_events_NSC_apn_ge10
        flag_star_detected flag_star_NSC_apn_gt0
        flag_star_NSC_apn_ge2 flag_star_NSC_apn_ge5 flag_star_NSC_apn_ge10;
run;

proc sort data=junc_flags;
   by junction_id;
proc sort data=pb_junc ;
   by junction_id;
run;

data pb_junc_support_flags;
  merge junc_flags (in=in1) pb_junc (in=in2);
  by junction_id;
  if in1 and in2;
run;

/* collapse flags to transcript level */

proc sort data=pb_junc_support_flags;
  by pacbio_id;
proc means data=pb_junc_support_flags noprint;
  by pacbio_id;
  var flag_events_detected flag_events_NSC_apn_gt0 flag_events_NSC_apn_ge2
      flag_events_NSC_apn_ge5 flag_events_NSC_apn_ge10
      flag_star_detected flag_star_NSC_apn_gt0 flag_star_NSC_apn_ge2
      flag_star_NSC_apn_ge5 flag_star_NSC_apn_ge10;
  output out=pb_junc_support_by_xs max=;
run;

/* Merge with PacBio transcripts */

proc sort data=pb_junc_support_by_xs;
  by pacbio_id;
proc sort data=pb_novel_flag_novel_junc;
  by pacbio_id;
run;

data pb_novel_xs_junc_support;
   merge pb_novel_flag_novel_junc (in=in1) pb_junc_support_by_xs (in=in2);
   by pacbio_id;
   if not in2 then do;
         flag_events_detected=0;
         flag_events_NSC_apn_gt0=0;
         flag_events_NSC_apn_ge2=0;
         flag_events_NSC_apn_ge5=0;
         flag_events_NSC_apn_ge10=0;
         flag_star_detected=0;
         flag_star_NSC_apn_gt0=0;
         flag_star_NSC_apn_ge2=0;
         flag_star_NSC_apn_ge5=0;
         flag_star_NSC_apn_ge10=0;
         end;
   if in1 then output;
run;

/* Count transcripts */

proc freq data=pb_novel_xs_junc_support noprint;
   tables flag_multigene*flag_has_junction*flag_has_novel_junction*
          flag_events_detected*flag_star_detected / out=pb_novel_detect;
   tables flag_multigene*flag_has_junction*flag_has_novel_junction*
          flag_events_NSC_apn_gt0*flag_star_NSC_apn_gt0 / out=pb_novel_apn0;
   tables flag_multigene*flag_has_junction*flag_has_novel_junction*
          flag_events_NSC_apn_ge2*flag_star_NSC_apn_ge2 / out=pb_novel_apn2;
   tables flag_multigene*flag_has_junction*flag_has_novel_junction*
          flag_events_NSC_apn_ge5*flag_star_NSC_apn_ge5 / out=pb_novel_apn5;
   tables flag_multigene*flag_has_junction*flag_has_novel_junction*
          flag_events_NSC_apn_ge10*flag_star_NSC_apn_ge10 / out=pb_novel_apn10;
run;

/* Detection */

proc print data=pb_novel_detect;
run;

/*
                           flag_has_                      flag_
   flag_      flag_has_      novel_     flag_events_      star_
 multigene     junction     junction      detected      detected    COUNT

     0            0            0              0             0         172
     0            1            0              0             0        1149

     0            1            1              0             0        1694
     0            1            1              0             1         340
     0            1            1              1             0          57
     0            1            1              1             1         145

     1            0            0              0             0         604
     1            1            0              0             0         915

     1            1            1              0             0        1002
     1            1            1              0             1         215
     1            1            1              1             0          30
     1            1            1              1             1         101

*/


/* APN>0 support */

proc print data=pb_novel_apn0;
run;

/*

                          flag_has_                    flag_star_
  flag_      flag_has_      novel_     flag_events_     NSC_apn_
multigene     junction     junction     NSC_apn_gt0        gt0       COUNT

    0            0            0              0              0          172
    0            1            0              0              0         1149

    0            1            1              0              0         1964
    0            1            1              0              1          155
    0            1            1              1              0           47
    0            1            1              1              1           70

    1            0            0              0              0          604
    1            1            0              0              0          915
    1            1            1              0              0         1137
    1            1            1              0              1          117
    1            1            1              1              0           35
    1            1            1              1              1           59


*/


/* APN>=2 support */

proc print data=pb_novel_apn2;
run;


/*

                          flag_has_                    flag_star_
  flag_      flag_has_      novel_     flag_events_     NSC_apn_
multigene     junction     junction     NSC_apn_ge2        ge2       COUNT

    0            0            0              0              0          172
    0            1            0              0              0         1149
    0            1            1              0              0         2189

    0            1            1              0              1            9
    0            1            1              1              0           38

    1            0            0              0              0          604
    1            1            0              0              0          915

    1            1            1              0              0         1304
    1            1            1              0              1            3
    1            1            1              1              0           41



*/

/* APN>=5 support */

proc print data=pb_novel_apn5;
run;


/*

                          flag_has_                    flag_star_
  flag_      flag_has_      novel_     flag_events_     NSC_apn_
multigene     junction     junction     NSC_apn_ge5        ge5       COUNT

    0            0            0              0              0          172
    0            1            0              0              0         1149

    0            1            1              0              0         2217
    0            1            1              0              1            3
    0            1            1              1              0           16

    1            0            0              0              0          604
    1            1            0              0              0          915

    1            1            1              0              0         1326
    1            1            1              1              0           22


*/

/* APN>=10 support */

proc print data=pb_novel_apn10;
run;


/*

                          flag_has_                    flag_star_
  flag_      flag_has_      novel_     flag_events_     NSC_apn_
multigene     junction     junction    NSC_apn_ge10       ge10       COUNT

    0            0            0              0              0          172
    0            1            0              0              0         1149

    0            1            1              0              0         2225
    0            1            1              1              0           11
    1            0            0              0              0          604

    1            1            0              0              0          915

    1            1            1              0              0         1335
    1            1            1              1              0           13

*/
