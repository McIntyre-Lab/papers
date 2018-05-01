/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Count number of matching exons per Refseq transcript */

data sim_exons2mm10; 
   set eventloc.simulated_exons_to_mm10;
run;

proc sort data=sim_exons2mm10;
   by simulation test sim_gene_id transcript_id;
proc freq data=sim_exons2mm10 noprint;
   by simulation test sim_gene_id;
   tables transcript_id / out=mm10_ex_per_sim_gene;
run;

data sim_exon_info;
   set sim_exons2mm10;
   keep simulation test sim_gene_id num_exons_per_sim_gene;
run;

proc sort data=mm10_ex_per_sim_gene ;
   by simulation test sim_gene_id ;
proc sort data=sim_exon_info nodup;
   by simulation test sim_gene_id ;
run;

data mm10_ex_per_sim_gene2;
  merge sim_exon_info (in=in1) mm10_ex_per_sim_gene (in=in2);
  by simulation test sim_gene_id ;
  if in1 and in2;
  drop PERCENT;
  rename COUNT=num_matching_exons;
run;

data mm10_exon_info;
   set sim_exons2mm10;
   keep transcript_id num_exons_per_xscript;
run;

proc sort data=mm10_ex_per_sim_gene2;
   by transcript_id ;
proc sort data=mm10_exon_info nodup;
   by transcript_id;
run;

data mm10_ex_per_sim_gene3;
  merge mm10_exon_info (in=in1) mm10_ex_per_sim_gene2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Flag if all exons for simulated genes match, and match to all exons for each possible Refseq transcript */

data flag_match;
  set mm10_ex_per_sim_gene3;
  if num_matching_exons = num_exons_per_sim_gene then flag_all_sim_match = 1;
  else flag_all_sim_match = 0;
  if num_matching_exons = num_exons_per_xscript then flag_all_mm10_match = 1;
  else flag_all_mm10_match = 0;
run;

proc freq data=flag_match;
  tables flag_all_sim_match*flag_all_mm10_match ;
run;

/*
flag_all_sim_match
          flag_all_mm10_match

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |2227970 |   1310 |2229280
         |  86.01 |   0.05 |  86.06
         |  99.94 |   0.06 |
         |  86.54 |   8.29 |
---------+--------+--------+
       1 | 346585 |  14487 | 361072
         |  13.38 |   0.56 |  13.94
         |  95.99 |   4.01 |
         |  13.46 |  91.71 |
---------+--------+--------+
Total     2574555    15797  2590352
            99.39     0.61   100.00

*/
/* Calculate the proportion of RefSeq exons matching and proportion of simulated exons matching.
   Then, I want to bin matches into the following sets:
   (1) All Refseq exons match and all simulated exons match (complete match)
   (2) All Refseq exons match, only 1 simulated exon is missin (complete refseq match)
   (3) All simulated exons match, only 1 Refseq exon is missing (complete simulated match)
   (4) 1 sim and 1 refseq exon not matching (partial match)
   Else, no match and delete observation
   
   Then sort by proportion of RefSeq exons matching and flag as best match */

data prop_matching;
   set flag_match;
   prop_refseq_match=num_matching_exons/num_exons_per_xscript;
   prop_sim_match=num_matching_exons/num_exons_per_sim_gene;
run;

data bin_matches;
   set prop_matching;
   if prop_refseq_match=1 then flag_full_refseq_match=1; else flag_full_refseq_match=0;
   if prop_sim_match=1 then flag_full_sim_match=1; else flag_full_sim_match=0;
   if num_exons_per_sim_gene - num_matching_exons = 1 then flag_sim_one_missing=1;
   else flag_sim_one_missing=0;
   if num_exons_per_xscript - num_matching_exons = 1 then flag_refseq_one_missing=1;
   else flag_refseq_one_missing=0;
run;

proc freq data=bin_matches noprint;
   tables flag_full_refseq_match*flag_full_sim_match*flag_sim_one_missing*flag_refseq_one_missing
          / out=exon_match_counts;
proc print data=exon_match_counts;
run;

/*
 flag_full_
   refseq_     flag_full_     flag_sim_     flag_refseq_
    match       sim_match    one_missing     one_missing      COUNT
      1             1             0               0           14487
      1             0             1               0             943
      1             0             0               0             367
      0             1             0               1           10536
      0             1             0               0          336049
      0             0             1               1           32198
      0             0             0               0         1523659
      0             0             0               1           13272
      0             0             1               0          658841
*/

data keep_matches;
  set bin_matches;
  if flag_full_refseq_match=1 and flag_full_sim_match=1 then output;
  if flag_full_refseq_match=1 and flag_sim_one_missing=1 then output;
  if flag_full_sim_match=1 and flag_refseq_one_missing=1 then output;
  if flag_sim_one_missing=1 and flag_refseq_one_missing=1 then output;
run;

proc sort data=keep_matches;
  by simulation test sim_gene_id descending prop_refseq_match;
run;

data flag_best_match;
   set keep_matches;
   by simulation test sim_gene_id;
   if first.sim_gene_id then flag_best_refseq_hit=1;
   else do;
     if flag_full_refseq_match=1 and flag_full_sim_match=1 then flag_best_refseq_hit=1;
     else flag_best_refseq_hit=0;
   end;
run;

proc sort data=flag_best_match;
   by simulation test;
proc freq data=flag_best_match noprint;
  by simulation test;
  tables flag_full_refseq_match*flag_full_sim_match*flag_sim_one_missing*flag_refseq_one_missing*flag_best_refseq_hit
   / out=match_summary;
proc print data=match_summary;
run;

/* Make best hits permenant so I can merge the Event analysis results with this */

data event.simulation_exon2refseq_best_hits;
   set flaG_best_match;
run;

