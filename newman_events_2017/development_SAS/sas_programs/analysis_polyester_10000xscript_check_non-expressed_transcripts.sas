/* check reads from one transcript with 0% of features detected */


data xs_no_cov;
   set event.bin_xs_by_dtct_apn0_10k;
   where perc_features_dtct = 0;
run;

data xs1000;
  set event.polyester_xs_list_10k;
run;

proc sort data=xs_no_cov;
  by transcript_id;
proc sort data=xs1000;
  by transcript_id;
run;

data xs10000_no_cov;
  merge xs10000 (in=in1) xs_no_cov (in=in2);
  by transcript_id;
  if in1 and in2;
run; 

proc print data=xs10000_no_cov (keep=transcript_id);
run;



/* 12 transcripts missed at APN>0:

 transcript_
 id

 NR_039558
 XR_379682
 XR_396084
 XR_404885
 XR_861692
 XR_864732
 XR_867906
 XR_874883
 XR_876445
 XR_877351
 XR_884124
 XR_887809

Check: number of reads
check where reads are aligning
Check coverage of fragments/fusions

Start with NR_039558: First 10 reads in sample 1 map to:

chr12	36816223
chr12	36816205
chr12	36816223
chr12	36816205
chr12	36816223
chr12	36816205
chr12	36816223
chr12	36816205
chr12	36816205
chr12	36816223

SAM flags are 81,161,97,145, all indicate that at least one read of the read pair mapped
to the exon

NR_039558 has one exon:
chr12 36816205	36816278	100628578:1

chr12	36816204	36816278	S54081_SI

Corresponding fragment/fusion is S54081_SI
No coverage??

100628578:1

SAMtools by default when generating mpileups, ignores "anomolous read pairs", for example, where only
one of the reads maps but the other doesn't

*/
