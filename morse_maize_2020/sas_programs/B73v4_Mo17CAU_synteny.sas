libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";


/*
are expressed genes syntenic between B73 and Mo17? 

    input = B73v4 vs Mo17 CAU synteny from Maggie Woodhouse
 
    YES = synfind matches synmap
    no = 
*/

/* import synteny list from maggie */
proc import datafile = "!MCLAB/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/B73v4.36_Mo17CAU_synfind_synmap_02amm.txt"
out = syn
dbms = tab replace ;
guessingrows = MAX ;
run;

proc sort data = syn nodups ;
by _all_ ;
run;

data syn2 ;
set syn ;
length match_synfind_synmap $16. ;
if Mo17_CAU_Synfind_default = '' and  Mo17_CAU_Synmap_megablast_0_001_ = '' then match_synfind_synmap = 'neither' ;
else if Mo17_CAU_Synfind_default =  Mo17_CAU_Synmap_megablast_0_001_ then match_synfind_synmap = 'match' ;

else if Mo17_CAU_Synfind_default = '' and Mo17_CAU_Synmap_megablast_0_001_ ne '' then  match_synfind_synmap = 'only synmap' ;
else if Mo17_CAU_Synmap_megablast_0_001_ = '' and Mo17_CAU_Synfind_default ne '' then  match_synfind_synmap = 'only synfind' ;
else if Mo17_CAU_Synfind_default ne Mo17_CAU_Synmap_megablast_0_001_ then  match_synfind_synmap = 'no match ' ;
else   match_synfind_synmap = 'oops' ;
run ;  /* 49373 obs */

data find_oops ;
set syn2 ;
where  match_synfind_synmap = 'oops' ;
run;  /* 0 in oops */


proc freq data = syn2 noprint;
tables  B73v4_gene_ID / out = cnts ;
run;  /* 49202 genes */

data ck ;
set cnts ;
where count ne 1 ;
run;  /* 169 genes duplicated 2-3 times */

data flag_dups ;
set ck ;
flag_dup = 1 ;
keep B73v4_gene_ID flag_dup ;
run;

proc sort data = syn2 ;
by B73v4_gene_ID ;
proc sort data = flag_dups ;
by B73v4_gene_ID ;

data add_flag ;
merge flag_dups syn2  ;
by B73v4_gene_ID ;
run ;


data find_dups ;
set add_flag ;
where flag_dup = 1 ;
run;
proc freq data = find_dups ;
tables match_synfind_synmap ;
run;  /* all are match or no match */

proc sort data = add_flag ;
by B73v4_gene_ID match_synfind_synmap ;
run ;

data synteny ;
set add_flag ;
by B73v4_gene_ID ;
if first.B73v4_gene_ID ;
run;  /* 49202 obs ! */



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

num_genotypes = (on_B73 + on_C123 + on_Hp301 + on_Mo17 + on_NC338) ;

keep PBgeneID num_genotypes ;
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


proc sort data = subset_onCalls ;
by B73v4_gene_ID ;
proc sort data = synteny ;
by B73v4_gene_ID ;
run ;

data synteny2 missing;
merge subset_onCalls (in=in1) synteny (in=in2) ;
by B73v4_gene_ID ;
if in1 and in2 then output synteny2 ;
else if in1 and not in2 then output missing ;
run;

data synteny_ck ;
set synteny2 ;
length syntenic $4. ;
if match_synfind_synmap = 'match' then syntenic = 'yes' ;
else if match_synfind_synmap = 'only synfind' or match_synfind_synmap = 'only synmap' then syntenic = 'yes' ;

else if match_synfind_synmap  = 'neither' then syntenic = 'no' ;
else if match_synfind_synmap  = 'no match' then syntenic = 'no' ;
else syntenic = 'oops' ;
run;

data find ;
set synteny_ck ;
where syntenic = 'oops';
run;

proc freq data = synteny_ck ;
tables match_synfind_synmap syntenic num_genotypes ;
tables match_synfind_synmap * num_genotypes ;
run ;
/*
match_
synfind_
synmap          Frequency     Percen
------------------------------------
match              10546       83.22
neither              660        5.21
no match             243        1.92
only synfind        1071        8.45
only synmap          152        1.20



  syntenic    Frequency     Percent
  ----------------------------------
  no               903        7.13
  yes            11769       92.87

*/

data PacBio_B73v4_Mo17Cau_synteny ;
retain B73v4_gene_ID syntenic num_genotypes ;
set synteny_ck ;
rename B73v4_gene_ID = PacBio_ZMgn ;
rename num_genotypes = num_genotype_expressed ;
keep B73v4_gene_ID num_genotypes syntenic ;
run ;

proc  freq data = PacBio_B73v4_Mo17Cau_synteny ;
tables syntenic num_genotype_expressed syntenic * num_genotype_expressed ;
run;
/*
syntenic     num_genotype_expressed

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|       2|       3|       4|       5|  Total
---------+--------+--------+--------+--------+--------+--------+
no       |      1 |     26 |     52 |     50 |     61 |    713 |    903
         |   0.01 |   0.21 |   0.41 |   0.39 |   0.48 |   5.63 |   7.13
         |   0.11 |   2.88 |   5.76 |   5.54 |   6.76 |  78.96 |
         |   8.33 |  45.61 |  41.60 |  28.25 |  18.37 |   5.96 |
---------+--------+--------+--------+--------+--------+--------+
yes      |     11 |     31 |     73 |    127 |    271 |  11256 |  11769
         |   0.09 |   0.24 |   0.58 |   1.00 |   2.14 |  88.83 |  92.87
         |   0.09 |   0.26 |   0.62 |   1.08 |   2.30 |  95.64 |
         |  91.67 |  54.39 |  58.40 |  71.75 |  81.63 |  94.04 |
---------+--------+--------+--------+--------+--------+--------+
Total          12       57      125      177      332    11969    12672
             0.09     0.45     0.99     1.40     2.62    94.45   100.00
*/


proc export data = PacBio_B73v4_Mo17Cau_synteny 
outfile = "!MCLAB/maize_ozone_FINAL/pacbio_paper/tables/PacBio_B73v4_Mo17Cau_synteny.tsv"
dbms = tab replace ;
run ;


