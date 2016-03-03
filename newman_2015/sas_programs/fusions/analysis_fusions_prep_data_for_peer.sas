/* Get a list of fusions for calculating PEER factors */

libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname con '/home/jrbnewman/concannon/sas_data';

/* We only want to do this for fusions that are on in all cell types */

data all_on_fusions;
  set con.fusions_on_gt_apn0;
  if flag_fusion_all_on0=1;
  keep fusion_id;
run;

/* Now get the normalized expression data for these fusions -- BUT only the 241 samples that with good coverage! */

data fusion_counts;
  set con.fusion_q3_norm_data_all;

/* mhalanobis distance and low coverage*/
if name =  '2009-PC-0221' then delete; *sample 75 cd8;
if name =  '2009-PC-0144' then delete; *sample 48 cd4;
if name =  '2009-PC-0236' then delete; *sample 80;
if name =  '2009-PC-0237' then delete; *sample 80;
if name =  '2009-PC-0235' then delete; *sample 80;
run;

proc sort data=fusion_counts;
   by fusion_id;
proc sort data=all_on_fusions;
   by fusion_id;
run;

data fusion_counts_all_on;
   merge fusion_counts (in=in1) all_on_fusions (in=in2);
   by fusion_id;
   if in1 and in2;
run;

/* Transpose and export for PEER */

proc sort data=fusion_counts_all_on;
   by fusion_id name;
run;

proc transpose data=fusion_counts_all_on
    out=fusion_counts_tpose;
   by fusion_id;
   id name;
   var log_q3_q3_apn;
run;

/* Export for PEER */

data peer_data;
  set fusion_counts_tpose;
  drop fusion_id _NAME_ ;
run;

proc export data=peer_data outfile='/home/jrbnewman/concannon/peer/expr_data_for_peer.csv' dbms=csv replace; putnames=no; run;


