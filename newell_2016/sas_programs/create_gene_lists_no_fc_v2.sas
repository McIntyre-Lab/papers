/* I am going to create gene lists without using the fold change that Michelle used. I will use fdr .05 and .10. 

*/

libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";

libname ribo "S:\McIntyre_Lab\arbeitman\arbeitman_ribotag\sas_data";


/* I'm going to leave the old model counts here for checks*/
proc freq data=ribo.results_by_fusion;
  tables flag_fdr_: ;
  run;

/*Contrasts: contrast 1 = ipMale-ipFemale
            contrast 2 = ipMale-inputMale
            contrast 3 = ipFemale-inputFemale
            contrast 4 = inputMale-inputFemale
            contrast 5 = male-female
            */


/*Results of proc freq:

FDR .20: IPmale-IPfemale        5,385 
         IPmale-InputMale       23,770
         IPfemale-InputFemale   23,695
         InputMale-InputFemale  2,797


FDR .10: IPmale-IPfemale        2,455
         IPmale-InputMale       23,740
         IPfemale-InputFemale   23,658
         InputMale-InputFemale  1,311


FDR .05: IPmale-IPfemale        1,048
         IPmale-InputMale       23,706
         IPfemale-InputFemale   23,616
         InputMale-InputFemale  675

*/

data logrpkm;
  set ribo.results_by_fusion;
  keep fusion_id ipmale_mean_logrpkm ipfemale_mean_logrpkm inputfemale_mean_logrpkm inputmale_mean_logrpkm;
  if flag_fusion_all_on0 = 1;
  run;

data logrpkm2;
  set logrpkm;
  ip_diff = ipmale_mean_logrpkm - ipfemale_mean_logrpkm;
  input_diff = inputmale_mean_logrpkm - inputfemale_mean_logrpkm;
  female_diff = ipfemale_mean_logrpkm - inputfemale_mean_logrpkm;
  male_diff = ipmale_mean_logrpkm - inputmale_mean_logrpkm;
  diff = (ipmale_mean_logrpkm - inputmale_mean_logrpkm) - (ipfemale_mean_logrpkm - inputfemale_mean_logrpkm);
  run;

data logrpkm3;
  set logrpkm2;
  if ip_diff > 0 then flag_male_TRAP_higher = 1; else flag_maleTRAP_higher =0;
  if input_diff > 0 then flag_inputmale_higher = 1; else flag_inputmale_higher = 0;
  if female_diff > 0 then flag_female_trap_higher = 1; else flag_female_trap_higher = 0;
  if male_diff > 0 then flag_male_trap_higher =1; else flag_male_trap_higher =0;
  if diff > 0 then flag_male_bias = 1; else flag_male_bias = 0;
  run;

proc freq data=logrpkm3;
  tables flag_: ;
  run;

/* Results of proc freq (logrpkm):

Higher in male TRAP than female TRAP:                   16,524
Higher in male TRAP than male input:                    3,886
Higher in female TRAP than female input:                2,402
Higher in male input than female input:                 14,432
Higher in male than female TRAP (adjusted for input):   15,829
*/

proc freq data=logrpkm3;
  tables flag_maleTRAP_higher * flag_male_enrich ;
  run;
proc freq data=logrpkm3;
  tables flag_maleTRAP_higher * flag_male_enrich ;
  run;


proc freq data=ribo.results_by_fusion;
  tables flag_fdr_p_contrast_3_05 * flag_fdr_p_contrast_2_05;
  run;





/* NEW MODEL */
/* Using the results from the updated proc glimmix, I made FDR flags, results here: 
*/

proc freq data = ribo.flag_fdr_contrast_by_fusion_2;
  tables flag_fdr_: ;
  run;

/*Results of proc freq (Total fusions analyzed: 23,814):

FDR .20: IPmale-IPfemale        2,037 
         IPmale-InputMale       11,942
         IPfemale-InputFemale   11,401
         InputMale-InputFemale  4,132
         Male-Female            811

FDR .10: IPmale-IPfemale        864
         IPmale-InputMale       8,777
         IPfemale-InputFemale   7,778
         InputMale-InputFemale  2,332
         Male-Female            238

FDR .05: IPmale-IPfemale        388
         IPmale-InputMale       6,240
         IPfemale-InputFemale   5,281
         InputMale-InputFemale  1,493
         Male-Female            75
*/


proc freq data=ribo.results_by_fusion_new_model;
  tables flag_fdr_p_contrast_3_05 * flag_fdr_p_contrast_2_05;
  run;

/* Results from proc freq. Total fusions: 51,506 

Neither different from input:  15,196
Females different from input:  2,378
Males different from input:    3,337
Both different from input:     2,903

*/

/*Cross tab female and male ip enriched with male enriched, FDR 0.05*/
proc freq data=ribo.results_by_fusion_new_model;
  tables flag_fdr_p_contrast_5_05 * flag_fdr_p_contrast_2_05 * flag_fdr_p_contrast_3_05;
  run;

/*femaleip-femaleinput * maleip-maleinput controlling for male enriched

When male enriched = 1 (411)
    Males different from input:     41
    Females different from input:   17
    Neither different from input:   6
    Both different from input:      11

When male enriched = 0 (23,739)
    Males different from input:     2,337
    Females different from input:   3,3320
    Neither different from input:   15,190
    Both different from input:      2,892

*/


/*Make flags needed for counts */
data results_by_fusion;
  set ribo.results_by_fusion_new_model;
    diff = (ipmale_mean_logrpkm - inputmale_mean_logrpkm) - (ipfemale_mean_logrpkm - inputfemale_mean_logrpkm);
    female_diff = ipfemale_mean_logrpkm - inputfemale_mean_logrpkm;
    male_diff = ipmale_mean_logrpkm - inputmale_mean_logrpkm;
    input_sex_diff = inputmale_mean_logrpkm - inputfemale_mean_logrpkm;
    if diff > 0 then flag_male_bias = 1 ; else flag_male_bias = 0;
    if female_diff > 0 then flag_female_trap_higher = 1; else flag_female_trap_higher = 0;
    if male_diff > 0 then flag_male_trap_higher =1; else flag_male_trap_higher =0;
    if input_sex_diff > 0 then flag_input_male_bias = 1; else flag_input_male_bias =0;
    if input_sex_diff < 0 then flag_input_female_bias = 1; else flag_input_female_bias = 0;
   run;



data results_by_fusion_2;
  set results_By_fusion;
  if flag_female_Trap_higher = 1 and flag_Fdr_p_contrast_3_05 =1 then flag_female_trap_sig_5 =1; else flag_female_trap_sig_5=0;
  if flag_male_Trap_higher = 1 and flag_Fdr_p_contrast_2_05 =1 then flag_male_trap_sig_5 =1; else flag_male_trap_sig_5=0;
    if flag_female_Trap_higher = 1 and flag_Fdr_p_contrast_3_10 =1 then flag_female_trap_sig_10 =1; else flag_female_trap_sig_10=0;
  if flag_male_Trap_higher = 1 and flag_Fdr_p_contrast_2_10 =1 then flag_male_trap_sig_10 =1; else flag_male_trap_sig_10=0;

  if flag_female_Trap_higher = 1 and flag_Fdr_p_contrast_3_20 =1 then flag_female_trap_sig_20 =1; else flag_female_trap_sig_20=0;
  if flag_male_Trap_higher = 1 and flag_Fdr_p_contrast_2_20 =1 then flag_male_trap_sig_20 =1; else flag_male_trap_sig_20=0;
   run;



/*Look for input differences*/
proc freq data=results_by_fusion_2;
  tables flag_fdr_p_contrast_4_20 * flag_input_male_bias * flag_input_female_bias;
  run;

/*
When the input sex difference is not significant:
    Female input biased: 7,466
    Male input biased: 12,216

When the input sex difference is significant:
    Female input biased: 1,916
    Male input biased: 2,216
*/




*FDR05;
proc freq data=results_by_fusion;
  tables flag_female_trap_higher * flag_fdr_p_contrast_3_05;
  run; *260 are significantly different from input and higher than input;

proc freq data=results_by_fusion;
  tables flag_male_trap_higher * flag_fdr_p_contrast_2_05;
  run; *512 are bigger than input AND the difference is significant;

*FDR20;
proc freq data=results_by_fusion;
  tables flag_female_trap_higher * flag_fdr_p_contrast_3_20;
  run;  *501 are significantly different from input and higher than input;

proc freq data=results_by_fusion;
  tables flag_male_trap_higher * flag_fdr_p_contrast_2_20;
  run; *998 are significantly different from input and higher than input;





*New tables;
proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_05 * flag_female_trap_sig_5;
  run; * 4 fusions are female trap significant AND the sex difference is significant;
proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_10 * flag_female_trap_sig_10;
  run; * 11 fusions;
proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_20 * flag_female_trap_sig_20;
  run; * 33 fusions;

proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_05 * flag_male_trap_sig_5;
  run; * 10 fusions are male trap significant AND the sex different is significant;
proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_10 * flag_male_trap_sig_10;
  run; * 30 fusions;
proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_20 * flag_male_trap_sig_20;
  run; * 104 fusions;

/*Trap sex specific AND sex different is signif*/
proc freq data = results_by_fusion_2;
  tables flag_fdr_p_contrast_5_20 * flag_male_trap_sig_20 * flag_female_trap_sig_20;
  run;
/* Results from proc freq

When the sex diff is NOT significant (Total=23,003 fusions):
    Male only signif: 559
    Female only signif: 133
    Both signif: 335

When the sex diff is significant (Total=811 fusions):
    Male only signif: 87
    Female only signif: 16
    Both signif: 17
*/

data male_trap_enrich;
  set results_by_fusion_2;
  if flag_male_trap_sig_20 = 1 and flag_female_trap_sig_20 = 0 and flag_fdr_p_contrast_5_20=1;
  run; *87 fusions;
  proc sort data=male_trap_enrich nodupkey;
    by fbgn_cat;
    run;*87 genes;

data female_trap_enrich;
  set results_by_fusion_2;
  if flag_female_trap_sig_20 = 1 and flag_male_trap_sig_20 = 0 and flag_fdr_p_contrast_5_20=1;
  run; *16 fusions;
  proc sort data=female_trap_enrich nodupkey;
    by fbgn_cat;
    run; *16 genes;



*Only where the sex difference is significant;
proc freq data=results_By_fusion_2;
  where flag_fdr_p_contrast_5_05=1;
  tables flag_male_bias * flag_female_trap_sig_5 * flag_male_trap_sig_5;
  run;
proc freq data=results_By_fusion_2;
  where flag_fdr_p_contrast_5_10=1;
  tables flag_male_bias * flag_female_trap_sig_10 * flag_male_trap_sig_10;
  run;
proc freq data=results_By_fusion_2;
  where flag_fdr_p_contrast_5_20=1;
  tables flag_male_bias * flag_female_trap_sig_20 * flag_male_trap_sig_20;
  run;

/* Results when the sex diff is always signif:
FDR 0.05:
    Male bias=0 (total=14)
        Female trap sig, not male trap sig: 2 fusions
        Male trap sig, not female trap sig: 0
        Neither trap sig: 12
        Both trap sig: 0
    Male bias=1 (total=61)
        Male trap sig, not female trap sig: 8 fusions
        Female trap sig, not male trap sig: 0
        Neither trap sig: 51
        both trap sig: 2
FDR 0.10:
    Male bias =0 (total=46)
        Female trap sig, not male trap sig: 3
        Male trap sig, not female trap sig:0
        Neither trap sig: 42
        Both trap sig:1
    Male bias=1 (total=192)
        Female trap sig, not male trap sig: 0
        Male trap sig, not female trap sig: 22
        Neither trap sig: 163
        Both trap sig: 7
FDR 0.20
    Male bias=0 (total=174)
        Female trap sig, not male trap sig: 16
        Male trap sig, not female trap sig:0
        Neither trap sig: 155
        Both trap sig:3
    Male bias=1 (total=637)
        Female trap sig, not male trap sig: 0
        Male trap sig, not female trap sig: 87
        Neither trap sig: 536
        Both trap sig: 14      

*/


*Sex diff sig X male bias;
proc freq data=results_by_fusion_2;
  tables flag_fdr_p_contrast_5_05 * flag_male_bias;
  run;
proc freq data=results_by_fusion_2;
  tables flag_fdr_p_contrast_5_10 * flag_male_bias;
  run;
proc freq data=results_by_fusion_2;
  tables flag_fdr_p_contrast_5_20 * flag_male_bias;
  run;
/*
FDR 0.05 (23,814 fusions total)
    male biased and sex diff sig: 61
FDR 0.10
    male biasd and sex diff sig: 192
FDR 020
    male biased and sex diff sig: 637
*/




data trap_enriched;
  set results_by_fusion_2;
  if flag_female_trap_sig_20 = 1 and flag_male_trap_sig_20 = 1 and flag_fdr_p_contrast_1_20 = 0;
  run; *278 fusions;
proc sort data=trap_enriched nodupkey;
  by symbol_cat;
  run; *258 genes; 
  /* No FRU. The only fru fusion that is male and female trap sig, is also sex-diff-sig*/





/*Check for fru*/
data fru;
  set results_by_fusion_2;
  if symbol_cat = 'fru';
  run;
proc sort data =fru;
 by start;
 run;

proc freq data=ribo.on_calls_gt_apn0;
  tables flag_ipmale_on * flag_inputmale_on;
  tables flag_ipfemale_on * flag_inputfemale_on;
  run;

data check; 
  set ribo.on_calls_gt_apn0;
  if flag_ipmale_on=1 and flag_inputmale_on=0
  or flag_ipfemale_on=1 and flag_inputfemale_on=0;
  
  run;


