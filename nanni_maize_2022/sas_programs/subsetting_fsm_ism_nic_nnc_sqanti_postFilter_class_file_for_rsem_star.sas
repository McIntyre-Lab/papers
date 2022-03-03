libname pb "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*
(1) extract fsm and ism
  sort by zmtr, descending length, descending isoform (superPBID)
  pick longest superPBID for each zmtr
  if more that 1 with same length, pick larged superPBID number

    output:
        pb.fsm_ism_isoform (table with 19,251 PBIDs)
        pb.fsm_ism_isoform_zmtr (table with 19,251 PBIDs and ZMtrs)
        MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/fsm_ism_isoform.csv

(2) extract all nic and nnc
    
    output:
        pb.nic_nnc_isoform (table with 14,016 PBIDs)
        pb.nic_nnc_isoform_zmgn (table with 14,016 PBIDs and ZMgns)
        MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/nic_nnc_isoform.csv
   
*/ 

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_post_filter_b73/SQANTI_classification.txt"
out = class_postFilter_chained 
dbms = tab replace ;
guessingrows = max ;
run;   /* 48950 obs */


data class2_postFilter_chained ;
set class_postFilter_chained ;
if structural_category = "novel_not_in_catalog" or 
    structural_category = "incomplete-splice_match" or 
    structural_category = "full-splice_match" or 
    structural_category = "novel_in_catalog" 
    then flag_fsm_ism_nic_nnc = 1 ;
else flag_fsm_ism_nic_nnc = 0 ;
if structural_category = "incomplete-splice_match" or 
    structural_category = "full-splice_match" 
    then flag_fsm_ism = 1;
else flag_fsm_ism = 0 ;
if structural_category = "novel_not_in_catalog" or 
    structural_category = "novel_in_catalog" 
    then flag_nic_nnc = 1 ;
else flag_nic_nnc  = 0;
if find(associated_gene, "novelGene") > 0 then flag_novel_gene = 1;
    else flag_novel_gene = 0;
run;


proc freq data =  class2_postFilter_chained ;
tables  flag_fsm_ism_nic_nnc 
        flag_fsm_ism 
        flag_nic_nnc
        flag_novel_gene;
run ;
/*
   flag_fsm_ism_
         nic_nnc    Frequency
------------------------------
               0        5325
               1       43625

  flag_fsm_ism    Frequency
  ----------------------------
             0       19341
             1       29609

flag_nic_nnc    Frequency
-------------------------
           0       34934
           1       14016

 flag_novel_gene    Frequency
 -----------------------------
               0       46279
               1        2671
*/


data subset_postFilter_chained ;
set class2_postFilter_chained ;
if flag_fsm_ism = 1 ;
run ;

proc sort data = subset_postFilter_chained ;
by associated_transcript descending length descending isoform;
run;

data keep_fsm_ism ;
retain count ;
set subset_postFilter_chained ;
count + 1 ;
by associated_transcript descending length descending isoform;
if first.associated_transcript then count = 1;
run ;

data fsm_ism ;
set keep_fsm_ism ;
where count = 1 ;
keep isoform associated_transcript ;
rename isoform = superPBID 
    associated_transcript = ZMtr ;
run;  /* 19251 obs */

data 

proc sort data = fsm_ism nodups ;
by _all_ ;
run;  /* all uniq */

data pb.fsm_ism_isoform_zmtr ;
set fsm_ism ;
run ;

data pb.fsm_ism_isoform ;
set pb.fsm_ism_isoform_zmtr ;
keep isoform ;
run ;

proc export data = pb.fsm_ism_isoform 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/fsm_ism_isoform.csv"
dbms = csv replace ;
run ;


/* (2) NNC and NIC */
data not_fsm ;
set class2_postFilter_chained ;
if  flag_nic_nnc = 1  ;
run ;

proc freq data = not_fsm noprint;
tables isoform / out = cnts ;
run ;

data ck ;
set cnts ;
where count ne 1 ;
run;  /* uniq!! */

data pb.nic_nnc_isoform_zmgn ;
set not_fsm ;
keep isoform associated_gene ;
rename associated_gene = ZMgn ;
run ;

data pb.nic_nnc_isoform ;
set  pb.nic_nnc_isoform_zmgn ;
keep isoform ;
run;

proc export data = pb.nic_nnc_isoform 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/nic_nnc_isoform.csv"
dbms = csv replace ;
run ;


data catting ;
set pb.fsm_ism_isoform pb.nic_nnc_isoform ;
run;

proc sort data = catting nodups ;
by _all_ ;
run;  /* good uniq after setting together */


