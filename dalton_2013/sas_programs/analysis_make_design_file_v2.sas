libname fru '!MCLAB/Fru_network/sasdata';
/* 
 *I have an issue with the sample_ids. According to the spreadsheet there should
 *be a 2011-04-26 however I only have a file with 2011-05-03. So I am going to
 *replace 2011-04-26 in the desin file to match the 2011-05-03. I am assuming
 *that these are the same data because the dates are so close
 *
 * REVISIONS: 12/23/2011
 *              - Removed all other treatments, except for Adult heads since
 *                this is all of interest right now
 */

data fru.design_file;
    length trt $19. sample_id $19.;
    input trt $ sample_id $ ;
    datalines;
    AH_BerF 2011-05-03_5_TGACCA
    AH_BerF 2011-05-03_3_ACAGTG
    AH_BerF 2011-05-03_2_GCCAAT
    AH_BerM 2011-05-03_1_ATCACG
    AH_BerM 2011-05-03_6_CGATGT
    AH_BerM 2011-05-03_5_TTAGGC
    AH_CS 2011-05-03_1_ACTTGA
    AH_CS 2011-05-03_2_ACTTGA
    AH_CS 2011-05-03_3_ACTTGA
    AH_CS 2011-05-03_5_ACTTGA
    AH_dsxD 2011-05-03_1_CAGATC
    AH_dsxD 2011-05-03_6_ACTTGA
    AH_dsxD 2011-05-03_5_GATCAG
    AH_dsxNullF 2011-05-03_6_TAGCTT
    AH_dsxNullF 2011-05-03_2_ATCACG
    AH_dsxNullF 2011-05-03_5_CGATGT
    AH_dsxNullM 2011-05-03_6_TTAGGC
    AH_dsxNullM 2011-05-03_3_TGACCA
    AH_dsxNullM 2011-05-03_6_ACAGTG
    AH_FruP14_440 2011-05-03_1_GATCAG
    AH_FruP14_440 2011-05-03_2_GATCAG
    AH_FruP14_440 2011-05-03_3_GATCAG
    AH_FruP14_440 2011-05-03_8_GATCAG
    AH_FruW12_ChaM5 2011-05-03_1_TAGCTT
    AH_FruW12_ChaM5 2011-05-03_2_TAGCTT
    AH_FruW12_ChaM5 2011-05-03_3_TAGCTT
    AH_FruW12_ChaM5 2011-05-03_7_TAGCTT
    AH_Female_FruM(A) 2011-07-05_1_TGACCA
    AH_Female_FruM(A) 2011-07-05_2_ACAGTG
    AH_Female_FruM(A) 2011-07-05_3_GCCAAT
    AH_Male_FruM(A) 2011-07-05_1_ATCACG
    AH_Male_FruM(A) 2011-07-05_2_CGATGT
    AH_Male_FruM(A) 2011-07-05_3_TTAGGC
    AH_Female_FruM(B) 2011-07-05_2_TGACCA
    AH_Female_FruM(B) 2011-07-05_3_ACAGTG
    AH_Female_FruM(B) 2011-07-05_1_GCCAAT
    AH_Male_FruM(B) 2011-07-05_2_ATCACG
    AH_Male_FruM(B) 2011-07-05_3_CGATGT
    AH_Male_FruM(B) 2011-07-05_1_TTAGGC
    AH_Female_FruM(C) 2011-07-05_5_TAGCTT
    AH_Female_FruM(C) 2011-07-05_6_ATCACG
    AH_Female_FruM(C) 2011-07-05_5_CGATGT
    AH_Male_FruM(C) 2011-07-05_1_CAGATC
    AH_Male_FruM(C) 2011-07-05_2_ACTTGA
    AH_Male_FruM(C) 2011-07-05_3_GATCAG
    AH_CSFemale 2011-07-05_7_CAGATC
    AH_CSFemale 2011-07-05_7_ACTTGA
    AH_CSFemale 2011-07-05_8_GATCAG
    ;
    run;
