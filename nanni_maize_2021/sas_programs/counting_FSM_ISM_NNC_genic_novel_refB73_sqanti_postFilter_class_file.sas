


libname pb "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*  extract NIC, NNC, genic and novel pacbio isoforms (FSM and ISM extracted in different script)
    ** output DOES NOT contain PB isoforms for which there is a FSM or ISM 

input:  
    sqanti qc on post filter chained ref= B73

    keep NIC, NNC, genic and novel genes (FSM and ISM extracted in different script)

output:
    pb.pbID_nnc_nic_genic_novel
    MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/pbID_nnc_nic_genic_novel.csv
    

*/


proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_post_filter_b73/SQANTI_classification.txt"
out = class_postFilter_chained
dbms = tab replace ;
guessingrows = max ;
run;   /* 48950 obs */

title "postFilter_chained" ;
proc freq data = class_postFilter_chained ;
tables structural_category ;
run;

data class2_postFilter_chained ;
set class_postFilter_chained ;
if structural_category = "novel_not_in_catalog" or 
    structural_category = "incomplete-splice_match" or 
    structural_category = "full-splice_match" or 
    structural_category = "novel_in_catalog" or 
    structural_category = "genic"     
    then flag_fsm_ism_nic_nnc_gen = 1 ;
else flag_fsm_ism_nic_nnc_gen = 0 ;
if structural_category = "incomplete-splice_match" or 
    structural_category = "full-splice_match" 
    then flag_fsm_ism = 1 ;
else flag_fsm_ism = 0 ;
if structural_category = "novel_not_in_catalog" or 
    structural_category = "novel_in_catalog" or 
    structural_category = "genic"     
    then flag_nic_nnc_gen = 1 ;
else flag_nic_nnc_gen = 0 ;
if find(associated_gene, "novelGene") ge 1 
    then flag_novel_gene = 1;
else flag_novel_gene = 0 ;
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

proc freq data = class2_postFilter_chained ;
tables flag_nic_nnc_gen  flag_novel_gene flag_nic_nnc_gen * flag_novel_gene ;
run; /*
 flag_nic_nnc_gen
           flag_novel_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  29650 |   2671 |  32321
          |  60.57 |   5.46 |  66.03
          |  91.74 |   8.26 |
          |  64.07 | 100.00 |
 ---------+--------+--------+
        1 |  16629 |      0 |  16629
          |  33.97 |   0.00 |  33.97
          | 100.00 |   0.00 |
          |  35.93 |   0.00 |
 ---------+--------+--------+
 Total       46279     2671    48950
             94.54     5.46   100.00    */

/* look at PB isoforms per zm_geneID */
proc freq data = class2_postFilter_chained ;
tables associated_gene / out = cnts ;
run ;

ods graphics on ;
proc sgplot data = cnts ;
histogram count ;
yaxis ranges=(0-20 60-100) ;
xaxis ranges=(0-20 200-400) ;
run ;

proc univariate data = cnts ;
var count ;
run; /*
     Location                    Variability

 Mean     3.290314     Std Deviation            6.91853
 Median   2.000000     Variance                47.86603
 Mode     1.000000     Range                  351.00000
                       Interquartile Range      3.00000 */

/* look at PB isoforms per zm_transcriptID */
proc freq data = class2_postFilter_chained ;
tables associated_transcript / out = cnts ;
run ;

ods graphics on ;
proc sgplot data = cnts ;
histogram count ;
yaxis ranges=(0-20 60-100) ;
xaxis ranges=(0-20 200-400) ;
run ;

proc univariate data = cnts ;
var count ;
run; /*
         Location                    Variability

 Mean     2.542593     Std Deviation          139.43897
 Median   1.000000     Variance                   19443
 Mode     1.000000     Range                      19340
                       Interquartile Range            0
 */

/* look at zm_transcriptID per zm_genetID */
data zm_per_zm ;
set class2_postFilter_chained ;
keep associated_gene associated_transcript ;
run ;
proc sort data = zm_per_zm nodups ;
by _all_ ;
run;

proc freq data = zm_per_zm ;
tables associated_gene / out = cnts ;
run ;

ods graphics on ;
proc sgplot data = cnts ;
histogram count ;
yaxis ranges=(0-20 60-100) ;
xaxis ranges=(0-20 200-400) ;
run ;

proc univariate data = cnts ;
var count ;
run; /*    Location                    Variability

Mean     1.828258     Std Deviation            1.33872
Median   1.000000     Variance                 1.79216
Mode     1.000000     Range                   13.00000
                      Interquartile Range      1.00000  */



proc sort data = cnts ;
by descending count  ;
run;

data flag ;
set cnts ;
if count = 1 then flag_1 = 1 ; else flag_1 = 0 ;
if count = 2 then flag_2 = 1 ; else flag_2 = 0 ;
if count > 2 then flag_gr_2 = 1 ; else flag_gr_2 = 0 ;
if count > 10 then flag_gr_10 = 1 ; else flag_gr_10 = 0 ;
if count > 50 then flag_gr_50 = 1 ; else flag_gr_50 = 0 ;
if count > 100 then flag_gr_100 = 1 ; else flag_gr_100 = 0 ;
if count > 200 then flag_gr_200 = 1 ; else flag_gr_200 = 0 ;
keep associated_gene flag_1 flag_2 flag_gr_2 flag_gr_10 flag_gr_50 flag_gr_100 flag_gr_200 ;
run ;


proc freq data = flag ;
tables flag_1 flag_2 ;
run ; /*
 flag_1    Frequency
 --------------------
      0        7978
      1        6899



 flag_2    Frequency
 --------------------
      0       12246
      1        2631     */
%macro flag (cnts) ;

proc freq data = flag ;
tables flag_gr_&cnts.  ;
run ;
%mend ;

%flag (2) ;
%flag (10) ;
%flag (50) ;
%flag (100) ;
%flag (200) ;/*
 flag_gr_2    Frequency
 ----------------------
         0        9530
         1        5347
flag_gr_10    Frequency
-----------------------
         0       14175
         1         702
flag_gr_50    Frequency
------------------------
         0       14852
         1          25
flag_gr_100    Frequency
------------------------
          0       14867
          1          10
flag_gr_200    Frequency
------------------------
          0       14873
          1           4 */

/* genes with >1 structural category */
 
proc sort data = class2_postFilter_chained ;
by  associated_gene;
proc freq data = class2_postFilter_chained noprint;
by associated_gene ;
tables structural_category / out = category_ck ;
run;
proc freq data = category_ck noprint;
tables  associated_gene / out = geneID_category;
run;

proc freq data = geneID_category;
tables count ;
run; /*
 COUNT    Frequency     Percent
 ------------------------------
     1        8266       55.56
     2        3725       25.04
     3        2043       13.73
     4         776        5.22
     5          67        0.45  */

/* transcripts with >1 structural category */
 
proc sort data = class2_postFilter_chained ;
by  associated_transcript;
proc freq data = class2_postFilter_chained noprint;
by associated_transcript ;
tables structural_category / out = category_ck ;
run;
proc freq data = category_ck noprint;
tables  associated_transcript / out = transcriptID_category;
run;

proc freq data = transcriptID_category;
tables count ;
run; /* COUNT    Frequency     Percent
 ------------------------------
     1       16422       85.30
     2        2829       14.69
     7           1        0.01  */



/* subset to not fsm ism */
data not_fsm ;
set class2_postFilter_chained ;
if  flag_nic_nnc_gen = 1 or flag_novel_gene = 1 ;
rename associated_gene = zm_geneID 
        associated_transcript = zm_transcriptID ;
run ;

/* keep pbIDs for these and export for subsetting gtf file to only these */
data pbID_nnc_nic_genic_novel ;
set not_fsm ;
keep isoform ;
rename isoform = pbID;
run ;

proc freq data = pbID_nnc_nic_genic_novel noprint;
tables isoform / out = cnt_iso ;
run ;
data check ;
set cnt_iso ;
where count ne 1 ;
run;  /* uniq list of pb isoforms */






/* find ones that are FSM or ISM AND have NNC NIC transcripts as well -- ones with multiple structural categories*/
proc sort data = not_fsm ;
by zm_geneID ;
proc freq data = not_fsm noprint;
by zm_geneID ;
tables structural_category / out = category_ck ;
run;
proc freq data = category_ck noprint;
tables  zm_geneID / out = geneID_category;
run;

proc freq data = geneID_category;
tables count ;
run;/*  102 geneIDs with 3 structural categories 
        1327 with 2 structural categories
        6480 with only 1 structural category*/

/* merge in flag for which zm_geneIDs have fsm or ism */
data fsm_ism ;
set  Zmtr_Zmgn_FSM_ISM ;
flag_gene_fsm_ism = 1 ;
run ;

proc sort data = fsm_ism;
by zm_geneID ;
proc sort data = not_fsm_postFilter_chained ;
by zm_geneID ;
run;

data both one two;
merge not_fsm_postFilter_chained (in=in1) fsm_ism (in=in2)  ;
by zm_geneID ;
if in1 ;
run;

data not_fsm_postFilter_chained2 ;
set both ;
if flag_gene_fsm_ism ne 1 then flag_gene_fsm_ism = 0;
run;

proc freq data = not_fsm_postFilter_chained2 ;
tables flag_gene_fsm_ism ;
run ;
/*
 flag_gene_
    fsm_ism    Frequency     Percent
------------------------------------
          0        4442       20.52
          1       17201       79.48 
20% of isoforms are not associated with a gene and that has a fsm or ism match  */

data no_fsm_ism_isoform ;
set not_fsm ;
where flag_gene_fsm_ism = 0 ;
run;


/* for each isoform with NO FSM or ISM match, count # ZM_geneIDs with >1 transcript */
proc freq data = no_fsm_ism_isoform noprint ;
tables zm_geneID / out = NN_cnts ;
run;
data ck_nn_cnts ;
set NN_cnts ;
where count ne 1 ;
run ;  /* there are 342 zm_geneIDs with more than 1 zm_transcriptID   */




/* look at PB isoforms per zm_geneID */
proc freq data = not_fsm_postFilter_chained2 ;
tables zm_geneID / out = cnts ;
run ;

proc univariate data = cnts ;
var count ;
run ;

proc sort data = cnts ;
by descending count  ;
run;

%macro counts (cnts) ;

data flag ;
set cnts ;
if count > &cnts. then flag = 1 ;
keep Zm_geneID flag ;
run ;

proc sort data = flag ;
by zm_geneID ;
proc sort data = chain_class_files ;
by zm_geneID ;
run ;

data look_&cnts. ;
merge flag (in=in1) chain_class_files (in=in2) ;
by zm_geneID ;
run;

data look2_&cnts. ;
set look_&cnts. ;
where flag = 1 ;
run;

title "num genes with > &cnts. transcripts per gene";
proc freq data = look2_&cnts. ;
tables chrom ;
run;
title ;
%mend ;

%counts (200) ;
%counts (100) ;
%counts (50) ;








