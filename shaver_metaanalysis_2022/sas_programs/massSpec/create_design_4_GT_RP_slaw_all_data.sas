libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;

/*

design for Brie GT RP SLAW run -- RP POS

import 'design' created from split wide python code
import original design file

merge together

*/



/* import design files from data */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/text_data/meta_analysis_rp_pos_SLAW_output_0xBFF_PH_design.tsv"
out = dsgn_rp_pos
dbms = tab replace ;
guessingrows = max ;
run ;

data dsgn_rp_pos ;
set dsgn_rp_pos ;
drop var2 ;
run;


/* import orig design */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_batchfile_param02_targets.tsv"
out = orig 
dbms = tab replace ;
guessingrows = max ;
run ;


proc sort data = dsgn_rp_pos ;
by sampleID ;
proc sort data = orig ;
by sampleID ;
run;

/* design files for GT RP output from Brie slaw */
data dsgn_rp_pos2 ;
merge dsgn_rp_pos (in=in1) orig (in=in2) ;
by sampleID ;
if in1 and in2 then flag_both =1 ;
else if in1 then flag_data = 1 ;
else if in2 then flag_orig = 1;
run ;


data dsgn_rp_pos3 ;
retain flag_: ;
set dsgn_rp_pos2 ;
if batch = . then do ;
    batch = substr(sampleID, 2, 1) ;
    end ;
else batch = batch ;

if sampleType = '' and sample ne '' then sampleType = sample ;
else if sampleType = '' and find(sampleID, "pool_pd1074") ge 1 then sampleType = "pool_pd1074";
else if sampleType = '' and find(sampleID, "pool_mutant") ge 1 then sampleType = "pool_mutant" ;
else if sampleType = '' and find(sampleID, "pool_solera") ge 1 then sampleType = "pool_solera" ;
else if sampleType = '' and find(sampleID, "extraction") ge 1 then sampleType = "blank" ;
else if sampleType = "pool" then sampleType = "pool_batch" ;
else sampleType = sampleType ;

if newBatch = '' then newBatch = scan(sampleID, 1, '_') ;
else newBatch = newBatch ;

sampleType_batch = compress(sampleType||'_'||newBatch) ;

if sampleType ne "PD1074" then strainNoPD1074 = strain  ;

run;

data dsgn_rp_pos4 ;
set dsgn_rp_pos3 ;
if flag_orig = 1 then delete ;
drop flag_both flag_data flag_orig ;

if sampleType = "pool_pd1074" then strain_w_poolPD = "pool_pd1074" ;
else strain_w_poolPD = strain ;
run;

proc freq data = dsgn_rp_pos4 ;
tables sampleType * batch ;
run;

data dsgn_GT_RP_POS_slaw ;
set dsgn_rp_pos4 ;
run ;

data sweet16.dsgn_GT_RP_POS_slaw ;
set dsgn_GT_RP_POS_slaw ;
run ;

proc export data = sweet16.dsgn_GT_RP_POS_slaw
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_slaw.tsv"
dbms = tab replace ;
run ;

data pds_rp_pos ;
set sweet16.dsgn_GT_RP_POS_slaw ;
where sampleType = "PD1074" or sampleType = "pool_pd1074";
run ;

proc export data = pds_rp_pos
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_slaw_onlyPDs.tsv"
dbms = tab replace ;
run ;


data new_rp_pos ;
set sweet16.dsgn_GT_RP_POS_slaw ;
where strain_w_poolPD ne '' ;
run;

proc export data = new_rp_pos
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_slaw_strainWPD.tsv"
dbms = tab replace ;
run ;

/* add in pair info */
proc import datafile ="/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/design_gt_all_batches_w_pair.tsv"
out = pairs 
dbms = tab replace ;
run;

data pairs2 ;
set pairs ;
keep sampleID pair set;
run;

proc sort data = pairs2 ;
by sampleID ;
proc sort data = sweet16.dsgn_GT_RP_POS_slaw ;
by sampleID ;
run;

data dsgn_GT_RP_POS_pairs_slaw oops ;
merge sweet16.dsgn_GT_RP_POS_slaw (in=in1) pairs2 (in=in2) ;
by sampleID ;
if in1 then output dsgn_GT_RP_POS_pairs_slaw ;
else output oops ;
run;  /* 0 in oops */

data sweet16.dsgn_GT_RP_POS_pairs_slaw ;
set dsgn_GT_RP_POS_pairs_slaw ;
if pair = . then delete ;
drop strainNoPD1074 ;
run;

proc export data = sweet16.dsgn_GT_RP_POS_pairs_slaw
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_pairs_slaw.tsv"
dbms = tab replace ;
run ;


/*
%macro batching (batch) ;

data pds_rp_pos_b&batch ;
set pds_rp_pos ;
if batch = &batch ;
run ;

proc export data = pds_rp_pos_b&batch 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_slaw_onlyPDs_batch&batch..tsv"
dbms = tab replace ;
run ;

%mend ;

%batching (1) ;
%batching (2) ;
%batching (3) ;
%batching (4) ;
%batching (5) ;
%batching (6) ;


