libname pb "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*  extract 'all' the Zm transcripts from 'all' classification text files using B73 as reference

input:  
    sqanti qc on post filter chained ref= B73


keep only FSM and ISM 

output:
    pb.Zmtr_FSM_ISM_all_class_files.sasdbat

    MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/Zmtr_FSM_ISM_all_class_files.csv  
    MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/structural_category_cnts_from_classFiles.pdf
*/


ods pdf file = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/structural_category_cnts_from_classFiles.pdf" ;

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_post_filter_b73/SQANTI_classification.txt"
out = class_postFilter_chained 
dbms = tab replace ;
guessingrows = max ;
run;   /* 48950 obs */

title "postFilter_chained" ;
proc freq data = class_postFilter_chained ;
tables structural_category ;
run;

/* structural_category        Frequency
 -------------------------------------
 antisense                       279
 full-splice_match             19354
 fusion                           41
 genic                          2613
 genic_intron                    125
 incomplete-splice_match       10255
 intergenic                     2267
 novel_in_catalog              10684
 novel_not_in_catalog           3332
*/

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
if find(associated_gene, "novelGene") > 0 then flag_novel_gene = 1;
    else flag_novel_gene = 0;
run;


proc freq data = class2_postFilter_chained ;
tables associated_gene / out = cnt_genes ;
run;  /* 14,877 */

proc freq data = class2_postFilter_chained ;
tables associated_transcript / out = cnt_transcript ;
run;  /*19,252 */

proc freq data = class2_postFilter_chained ;
tables isoform / out = cnt_isoforms ;
run;  /*48,950 */

proc freq data =  class2_postFilter_chained ;
tables  flag_fsm_ism_nic_nnc 
        flag_fsm_ism 
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

 flag_novel_gene    Frequency
 -----------------------------
               0       46279
               1        2671
*/

data subset_postFilter_chained ;
set class2_postFilter_chained ;
if flag_fsm_ism = 1 ;
run ;

data fsm_ism_postFilter_chained ;
set subset_postFilter_chained ;
keep associated_gene associated_transcript ;
rename associated_gene = Zm_geneID ;
rename associated_transcript = Zm_transcriptID ;
run;

proc sort data = fsm_ism_postFilter_chained nodups ;
by _all_ ;
run;

title ;


data _null_ ;
file print notitles ;
put '# unique Zm transcripts across all classification files is ' nobs=;
stop ;
set  fsm_ism_postFilter_chained nobs=nobs;
run ;
title ;

ods pdf close ;
 

