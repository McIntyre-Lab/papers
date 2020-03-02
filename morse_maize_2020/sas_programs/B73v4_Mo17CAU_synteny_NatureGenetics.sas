libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";

/*
are expressed genes syntenic between B73 and Mo17? 

input = B73v4 vs Mo17 CAU synteny from author of CAU paper Extensive intraspecific gene order and gene structural variations between Mo17 and other maize genomes. Nature Genetics. 50,1289-1295.
 
*/

/* import synteny list from junpeng shi*/
data WORK.SYN_B73    ;
%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
infile "!MCLAB/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/Mo17_Table2_gene_list/1.B73_Syntenic_genes/19.B73-sytenic_gene.txt"
    delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat gene $19. ;
       format gene $19. ;
    input
                gene  $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

data syntenic_b73 ;
set syn_b73 ;
flag_syntenic_b73 = 1 ;
run;

proc sort data = syntenic_b73 nodups ;
by _all_ ;
run;

/* DE */

data DE_prep ;
set tappas.tappas_results_genes ;
if dea_result_b73 = "DE" then flag_de_b73 = 1 ;
    else flag_de_b73 = 0;
if dea_result_c123 = "DE" then flag_de_c123 = 1 ;
    else flag_de_c123 = 0;
if dea_result_hp301 = "DE" then flag_de_hp301 = 1 ;
    else flag_de_hp301 = 0;
if dea_result_mo17 = "DE" then flag_de_mo17 = 1 ;
    else flag_de_mo17 = 0;
if dea_result_nc338 = "DE" then flag_de_nc338 = 1 ;
    else flag_de_nc338 = 0;

num_genotype_DE = (flag_de_b73 + flag_de_c123 + flag_de_hp301 + flag_de_mo17 + flag_de_nc338) ;
keep gene  num_genotype_DE;
run;

data DE ;
set DE_prep ;
if find(gene, "-T1") ge 1 then B73v4_gene_ID = tranwrd(gene, "-T1", "") ;
else  B73v4_gene_ID = gene ;
drop gene ;
run ;   /* 12391 obs */


/* combine with syntenic  */
data syntenic_b73_2 ;
set syntenic_b73 ;
rename gene = B73v4_gene_ID;
run ;

proc sort data = DE ;
by B73v4_gene_ID ;
proc sort data = syntenic_b73_2 ;
by B73v4_gene_ID ;
run ;


data synteny2 ;
merge DE (in=in1) syntenic_b73_2 (in=in2) ;
by B73v4_gene_ID ;
if in1 and in2;
run ;

data synteny3 ;
set synteny2 ;
if flag_syntenic_b73 ne 1 then flag_syntenic_b73 = 0 ;
run ;

proc freq data = synteny3 ;
tables flag_syntenic_b73 * num_genotype_DE / out = cnts;
run ;
/*
  flag_syntenic_b73     num_genotype_DE

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|       2|       3|       4|       5|  Total
  ---------+--------+--------+--------+--------+--------+--------+
         0 |    399 |    110 |     75 |     47 |     40 |      0 |    671
           |   3.22 |   0.89 |   0.61 |   0.38 |   0.32 |   0.00 |   5.42
           |  59.46 |  16.39 |  11.18 |   7.00 |   5.96 |   0.00 |
           |   5.34 |   5.69 |   6.64 |   5.46 |   4.36 |   0.00 |
  ---------+--------+--------+--------+--------+--------+--------+
         1 |   7070 |   1824 |   1054 |    814 |    877 |     81 |  11720
           |  57.06 |  14.72 |   8.51 |   6.57 |   7.08 |   0.65 |  94.58
           |  60.32 |  15.56 |   8.99 |   6.95 |   7.48 |   0.69 |
           |  94.66 |  94.31 |  93.36 |  94.54 |  95.64 | 100.00 |
  ---------+--------+--------+--------+--------+--------+--------+
  Total        7469     1934     1129      861      917       81    12391
              60.28    15.61     9.11     6.95     7.40     0.65   100.00

*/

/* synteny for expressed gene (on) */


/*  list of genes to check for synteny - fsm, ism, nnc and nic subset */

data fsm_ism ;
set pacbio.fsm_ism_isoform_zmtr ;
ZMgn =compress(scan(associated_transcript, 1, '_'));
keep ZMgn isoform ;
run ;

proc sort data = fsm_ism nodups ;
by _all_ ;
run;

data nic_nnc ;
set pacbio.nic_nnc_isoform_zmgn ;
run ;

proc sort data = nic_nnc nodups ;
by _all_ ;
run;


data fsm_ism_nic_nnc ;
set fsm_ism nic_nnc ;
run;

proc sort data = fsm_ism_nic_nnc nodups ;
by _all_ ;
run ;  /* 33,267 obs */

proc freq data = fsm_ism_nic_nnc ;
tables ZMgn / out = zmgn_cnt ;
run;   /* 12,764 genes */

data subset ;
set fsm_ism_nic_nnc ;
var1 = scan(isoform, 1, '.') ;
var2 = scan(isoform, 2, '.') ;
var3 = scan(isoform, 3, '.') ;
PBgeneID = compress(var1||'.'||var2) ;
drop var1-var3 isoform ;
run ;

proc sort data = subset nodups ;
by _all_ ;
run;

/* some duplicates - diff pb genes to same zmgn likely from nic and nnc matches to zmtr.... */
proc freq data = subset ;
tables zmgn / out = chk ;
run ;

data find ;
set chk ;
where count ne 1 ;
run;

data find2 ;
set subset ;
where ZMgn = 'Zm00001d001923';
run;

/* gene onCalls */
data gene_onCalls ;
set pacbio.sub_geno_trt_gene_oncall_tpm0;
keep geneID
    flag_B73_Amb_on0 flag_B73_Ele_on0 
    flag_C123_Amb_on0 flag_C123_Ele_on0  
    flag_Hp301_Amb_on0 flag_Hp301_Ele_on0 
    flag_Mo17_Amb_on0 flag_Mo17_Ele_on0 
    flag_NC338_Amb_on0 flag_NC338_Ele_on0 ;
rename geneID = PBgeneID ;
run ;

data gene_oncalls2 ;
set gene_oncalls ;
flag_sum_B73 = flag_B73_Amb_on0 + flag_B73_Ele_on0 ;
flag_sum_C123 = flag_C123_Amb_on0 + flag_C123_Ele_on0 ;
flag_sum_Hp301 = flag_Hp301_Amb_on0 + flag_Hp301_Ele_on0 ;
flag_sum_Mo17 = flag_Mo17_Amb_on0 + flag_Mo17_Ele_on0 ;
flag_sum_NC338 = flag_NC338_Amb_on0 + flag_NC338_Ele_on0 ;
run;

data gene_oncalls3 ;
retain PBgeneID num_genotypes ;
set gene_oncalls2 ;
if flag_sum_B73 ge 1 then on_B73 = 1 ; else on_B73 = 0 ;
if flag_sum_C123 ge 1 then on_C123 = 1 ; else on_C123 = 0 ;
if flag_sum_Hp301 ge 1 then on_Hp301 = 1 ; else on_Hp301 = 0 ;
if flag_sum_Mo17 ge 1 then on_Mo17 = 1 ; else on_Mo17 = 0 ;
if flag_sum_NC338 ge 1 then on_NC338 = 1 ; else on_NC338 = 0 ;

num_genotypes_on = (on_B73 + on_C123 + on_Hp301 + on_Mo17 + on_NC338) ;

keep PBgeneID num_genotypes_on ;
run ;

proc freq data = gene_oncalls3 ;
tables num_genotypes ;
run;

proc sort data = gene_oncalls3 ;
by PBgeneID ;
proc sort data = subset ;
by PBgeneID ;
run ;

data sub_onCalls  ouch;
merge subset (in=in1) gene_oncalls3 (in=in2) ;
by PBgeneID ;
if in1 and in2 then output sub_onCalls ;
else output ouch ;
run;  /* 0 in ouch */

data subset_onCalls ;
set sub_onCalls ;
rename ZMgn = B73v4_gene_ID ;
run ;

proc sort data = syntenic_b73_2 ;
by B73v4_gene_ID ;
proc sort data = subset_onCalls ;
by B73v4_gene_ID ;
run;

data synteny_on2 ;
merge subset_onCalls (in=in1) syntenic_b73_2 (in=in2) ;
by B73v4_gene_ID ;
if in1 ;
run ;

data synteny_on3 ;
set synteny_on2 ;
if flag_syntenic_b73 ne 1 then flag_syntenic_b73 = 0 ;
run ;

proc freq data = synteny_on3 ;
tables flag_syntenic_b73 * num_genotypes_on / out = cnts;
run ; 
/*
flag_syntenic_b73     num_genotypes_on

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|       2|       3|       4|       5|  Total
 ---------+--------+--------+--------+--------+--------+--------+
        0 |     11 |     27 |     49 |     50 |     46 |    573 |    756
          |   0.09 |   0.21 |   0.38 |   0.39 |   0.36 |   4.45 |   5.87
          |   1.46 |   3.57 |   6.48 |   6.61 |   6.08 |  75.79 |
          |  50.00 |  42.19 |  36.57 |  26.46 |  13.57 |   4.73 |
 ---------+--------+--------+--------+--------+--------+--------+
        1 |     11 |     37 |     85 |    139 |    293 |  11550 |  12115
          |   0.09 |   0.29 |   0.66 |   1.08 |   2.28 |  89.74 |  94.13
          |   0.09 |   0.31 |   0.70 |   1.15 |   2.42 |  95.34 |
          |  50.00 |  57.81 |  63.43 |  73.54 |  86.43 |  95.27 |
 ---------+--------+--------+--------+--------+--------+--------+
 Total          22       64      134      189      339    12123    12871
              0.17     0.50     1.04     1.47     2.63    94.19   100.00
*/
