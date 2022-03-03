
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";



/* 
freqs   on rsem subset (fsm, ism nic and nnc)
        isforms and genes
        tpm and expected counts

*/


/* TPM > 5 */
proc contents data = pacbio.sub_geno_trt_gene_onCall_tpm5;run;
data on_amb_all_genotypes ;
set pacbio.sub_geno_trt_gene_onCall_tpm5;
where   flag_B73_Amb_on5_allreps = 1 and  
        flag_C123_Amb_on5_allreps = 1 and 
        flag_Hp301_Amb_on5_allreps = 1 and 
        flag_Mo17_Amb_on5_allreps = 1 and 
        flag_NC338_Amb_on5_allreps = 1 ;
run ;
    /* 5243 obs out of 12,702 genes on in all amb reps in all 5 genotypes */
data on_ele_all_genotypes ;
set pacbio.sub_geno_trt_gene_onCall_tpm5;
where   flag_B73_ele_on5_allreps = 1 and  
        flag_C123_ele_on5_allreps = 1 and 
        flag_Hp301_ele_on5_allreps = 1 and 
        flag_Mo17_ele_on5_allreps = 1 and 
        flag_NC338_ele_on5_allreps = 1 ;
run ;
    /* 6379 obs out of 12,702 genes on in all amb reps in all 5 genotypes  */


/* TPM > 0 GENES */
proc contents data = pacbio.sub_geno_trt_gene_onCall_tpm0;run;
data on_amb_all_genotypes ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
where   flag_B73_Amb_on0_allreps = 1 and  
        flag_C123_Amb_on0_allreps = 1 and 
        flag_Hp301_Amb_on0_allreps = 1 and 
        flag_Mo17_Amb_on0_allreps = 1 and 
        flag_NC338_Amb_on0_allreps = 1 ;
run ;
    /* 10,195 obs out of 12,702 genes on in all amb reps in all 5 genotypes */
data on_ele_all_genotypes ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
where   flag_B73_ele_on0_allreps = 1 and  
        flag_C123_ele_on0_allreps = 1 and 
        flag_Hp301_ele_on0_allreps = 1 and 
        flag_Mo17_ele_on0_allreps = 1 and 
        flag_NC338_ele_on0_allreps = 1 ;
run ;
    /* 10,732 obs out of 12,702 genes on in all ele reps in all 5 genotypes  */

proc contents data = pacbio.sub_geno_trt_gene_onCall_tpm0;run;
data on_amb_genotypes ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
where   flag_B73_Amb_on0 = 1 and  
        flag_C123_Amb_on0 = 1 and 
        flag_Hp301_Amb_on0 = 1 and 
        flag_Mo17_Amb_on0 = 1 and 
        flag_NC338_Amb_on0 = 1 ;
run ;
    /* 11,637 obs out of 12,702 on in 50% of amb reps for all 5 genotypes */
data on_ele_all_genotypes ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
where   flag_B73_ele_on0 = 1 and  
        flag_C123_ele_on0 = 1 and 
        flag_Hp301_ele_on0 = 1 and 
        flag_Mo17_ele_on0 = 1 and 
        flag_NC338_ele_on0 = 1 ;
run ;
    /* 11,852 obs out of 12,702 on in 50% of ele reps for all 5 genotypes  */
data on_all_genotypes_gene ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
where   flag_B73_ele_on0 = 1 and  
        flag_C123_ele_on0 = 1 and 
        flag_Hp301_ele_on0 = 1 and 
        flag_Mo17_ele_on0 = 1 and 
        flag_NC338_ele_on0 = 1 and 
        flag_B73_amb_on0 = 1 and  
        flag_C123_amb_on0 = 1 and 
        flag_Hp301_amb_on0 = 1 and 
        flag_Mo17_amb_on0 = 1 and 
        flag_NC338_amb_on0 = 1
;
run ;      /* 11,550 obs out of 12,702 genes on in 50% of ele and amb reps for all 5 genotypes  */





/* TPM > 0 TRANSCRIPTS */
proc contents data = pacbio.sub_geno_trt_isoform_onCall_tpm0;run;
data on_amb_all_genotypes ;
set pacbio.sub_geno_trt_isoform_onCall_tpm0;
where   flag_B73_Amb_on0_allreps = 1 and  
        flag_C123_Amb_on0_allreps = 1 and 
        flag_Hp301_Amb_on0_allreps = 1 and 
        flag_Mo17_Amb_on0_allreps = 1 and 
        flag_NC338_Amb_on0_allreps = 1 ;
run ;      /* 10,752 obs out of 33,267 on in all amb reps in all 5 genotypes */

data on_ele_all_genotypes ;
set pacbio.sub_geno_trt_isoform_onCall_tpm0;
where   flag_B73_ele_on0_allreps = 1 and  
        flag_C123_ele_on0_allreps = 1 and 
        flag_Hp301_ele_on0_allreps = 1 and 
        flag_Mo17_ele_on0_allreps = 1 and 
        flag_NC338_ele_on0_allreps = 1 ;
run ;      /* 11,741 obs out of 33,267 on in all ele reps in all 5 genotypes  */

proc contents data = pacbio.sub_geno_trt_gene_onCall_tpm0;run;
data on_amb_genotypes ;
set pacbio.sub_geno_trt_isoform_onCall_tpm0;
where   flag_B73_Amb_on0 = 1 and  
        flag_C123_Amb_on0 = 1 and 
        flag_Hp301_Amb_on0 = 1 and 
        flag_Mo17_Amb_on0 = 1 and 
        flag_NC338_Amb_on0 = 1 ;
run ;      /* 18,532 obs out of 33,267 on in 50% of amb reps for all 5 genotypes */

data on_ele_all_genotypes ;
set pacbio.sub_geno_trt_isoform_onCall_tpm0;
where   flag_B73_ele_on0 = 1 and  
        flag_C123_ele_on0 = 1 and 
        flag_Hp301_ele_on0 = 1 and 
        flag_Mo17_ele_on0 = 1 and 
        flag_NC338_ele_on0 = 1 ;
run ;      /* 19,657 obs out of 33,267 on in 50% of ele reps for all 5 genotypes  */



data on_all_genotypes ;
set pacbio.sub_geno_trt_isoform_onCall_tpm0;
where   flag_B73_ele_on0 = 1 and  
        flag_C123_ele_on0 = 1 and 
        flag_Hp301_ele_on0 = 1 and 
        flag_Mo17_ele_on0 = 1 and 
        flag_NC338_ele_on0 = 1 and 
   flag_B73_amb_on0 = 1 and  
        flag_C123_amb_on0 = 1 and 
        flag_Hp301_amb_on0 = 1 and 
        flag_Mo17_amb_on0 = 1 and 
        flag_NC338_amb_on0 = 1
;
run ;      /* 17,664 obs out of 33,267 transcripts on in 50% of ele and amb reps for all 5 genotypes  */








proc contents data=pacbio.sub_geno_trt_isoform_onCall_tpm5; run ;

ods pdf file = "!MCLAB/maize_ozone_FINAL/2018/PacBio/text_data/freqs_rsem_subset_4_TPM_and_expected_counts.pdf" ;

%macro freqing (type, level) ;

title "&type:  genotype x trt,  TPM > &level. " ;
proc freq data=pacbio.sub_geno_trt_&type._onCall_tpm&level.;
tables flag_B73_Amb_on&level._allreps*flag_B73_ele_on&level._allreps ;
tables flag_C123_Amb_on&level._allreps*flag_C123_ele_on&level._allreps ;
tables flag_Hp301_Amb_on&level._allreps*flag_Hp301_ele_on&level._allreps ;
tables flag_Mo17_Amb_on&level._allreps*flag_Mo17_ele_on&level._allreps ;
tables flag_NC338_Amb_on&level._allreps*flag_NC338_ele_on&level._allreps ;

tables flag_B73_Amb_off&level._allreps*flag_B73_ele_off&level._allreps ;
tables flag_C123_Amb_off&level._allreps*flag_C123_ele_off&level._allreps ;
tables flag_Hp301_Amb_off&level._allreps*flag_Hp301_ele_off&level._allreps ;
tables flag_Mo17_Amb_off&level._allreps*flag_Mo17_ele_off&level._allreps ;
tables flag_NC338_Amb_off&level._allreps*flag_NC338_ele_off&level._allreps ;

run;

title "&type:  genotype x trt,  expected_count > &level. " ;
proc freq data=pacbio.sub_geno_trt_&type._onCall_Cnt&level.;
tables flag_B73_Amb_on&level._allreps*flag_B73_ele_on&level._allreps ;
tables flag_C123_Amb_on&level._allreps*flag_C123_ele_on&level._allreps ;
tables flag_Hp301_Amb_on&level._allreps*flag_Hp301_ele_on&level._allreps ;
tables flag_Mo17_Amb_on&level._allreps*flag_Mo17_ele_on&level._allreps ;
tables flag_NC338_Amb_on&level._allreps*flag_NC338_ele_on&level._allreps ;

tables flag_B73_Amb_off&level._allreps*flag_B73_ele_off&level._allreps ;
tables flag_C123_Amb_off&level._allreps*flag_C123_ele_off&level._allreps ;
tables flag_Hp301_Amb_off&level._allreps*flag_Hp301_ele_off&level._allreps ;
tables flag_Mo17_Amb_off&level._allreps*flag_Mo17_ele_off&level._allreps ;
tables flag_NC338_Amb_off&level._allreps*flag_NC338_ele_off&level._allreps ;
run;

title "&type:  genotype x trt, 50% of reps  TPM > &level. " ;
proc freq data=pacbio.sub_geno_trt_&type._onCall_tpm&level.;
tables flag_B73_Amb_on&level.*flag_B73_ele_on&level. ;
tables flag_C123_Amb_on&level.*flag_C123_ele_on&level.;
tables flag_Hp301_Amb_on&level.*flag_Hp301_ele_on&level.;
tables flag_Mo17_Amb_on&level.*flag_Mo17_ele_on&level.;
tables flag_NC338_Amb_on&level.*flag_NC338_ele_on&level.;

run;

title "&type:  genotype x trt, 50% of reps  expected_count > &level. " ;
proc freq data=pacbio.sub_geno_trt_&type._onCall_Cnt&level.;
tables flag_B73_Amb_on&level.*flag_B73_ele_on&level. ;
tables flag_C123_Amb_on&level.*flag_C123_ele_on&level.;
tables flag_Hp301_Amb_on&level.*flag_Hp301_ele_on&level.;
tables flag_Mo17_Amb_on&level.*flag_Mo17_ele_on&level.;
tables flag_NC338_Amb_on&level.*flag_NC338_ele_on&level.;
run;


title "&type:  genotype x trt,  TPM > &level. " ;
proc freq data = pacbio.sub_geno_trt_&type._onCall_tpm&level.;
tables flag_: ;
run;
    


title "&type:  genotype x trt,  expected_count > &level. " ;
proc freq data = pacbio.sub_geno_trt_&type._onCall_Cnt&level.;
tables flag_: ;
run;

%mend ;

%freqing (isoform, 0) ;
%freqing (isoform, 5) ;
%freqing (gene, 0) ;
%freqing (gene, 5) ;

ods pdf close ;
