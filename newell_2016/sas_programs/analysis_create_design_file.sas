/* Create design file for ribotag with sample ID info. trt ex = IPmale
*/


libname ribo "/home/fnew/mclab/arbeitman/arbeitman_ribotag/sas_data" ;


data ribo.design_file;
    length trt $19. type $19. sex $19. sample_id $29.;
    input trt $ type $ sex $ sample_id $;
    datalines;
InputFemale Female input Input_Female1_TTAGGC_L001_R1
InputFemale Female input Input_Female2_CAGATC_L001_R1
InputFemale Female input Input_Female3_GGCTAC_L002_R1
InputFemale Female input Input_Female4_CCGTCC_L002_R1
InputFemale Female input Input_Female5_GATCAG_L002_R1
InputMale Male input Input_Male1_ATCACG_L001_R1
InputMale Male input Input_Male2_ACAGTG_L001_R1
InputMale Male input Input_Male3_AGTCAA_L001_R1
InputMale Male input Input_Male4_GAGTGG_L001_R1
InputMale Male input Input_Male5_ATCACG_L002_R1
InputMale Male input Input_Male6_ACAGTG_L002_R1
IPfemale Female ip IP_Female1_TGACCA_L001_R1
IPfemale Female ip IP_Female2_ACTTGA_L001_R1
IPfemale Female ip IP_Female3_CTTGTA_L002_R1
IPfemale Female ip IP_Female4_GTCCGC_L002_R1
IPfemale Female ip IP_Female5_TAGCTT_L002_R1
IPmale Male ip IP_Male1_CGATGT_L001_R1
IPmale Male ip IP_Male2_GCCAAT_L001_R1
IPmale Male ip IP_Male3_AGTTCC_L001_R1
IPmale Male ip IP_Male4_ACTGAT_L002_R1
IPmale Male ip IP_Male5_CGATGT_L002_R1
IPmale Male ip IP_Male6_GCCAAT_L002_R1
;
run;
