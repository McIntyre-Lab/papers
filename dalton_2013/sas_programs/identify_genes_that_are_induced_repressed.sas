/********************************************************************************
* Using Venn Diagrams, I have identified several genes that are induced and
* repressed. Since we are summing fusions to give gene level results, it is
* likly a gene has one fusion that is induced and one fusion that is repressed.
* This would make since if one of the fusions is an Alternative exon. The goal
* of this script is to see if this is true.
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';


/* Create a list of genes that were induced and repressed */
%macro merge_ind_rep(sex,fru);

    data &sex._fruA_ind &sex._fruB_ind &sex._fruC_ind;
        set FRU.&sex._induced_union;
        if flag_a_ind = 1 then output &sex._fruA_ind;
        if flag_b_ind = 1 then output &sex._fruB_ind;
        if flag_c_ind = 1 then output &sex._fruC_ind;
        run;

    data &sex._fruA_rep &sex._fruB_rep &sex._fruC_rep;
        set FRU.&sex._repressed_union;
        if flag_a_rep = 1 then output &sex._fruA_rep;
        if flag_b_rep = 1 then output &sex._fruB_rep;
        if flag_c_rep = 1 then output &sex._fruC_rep;
        run;

    proc sort data = &sex._&fru._ind;
        by primary_fbgn;
        run;
    proc sort data = &sex._&fru._rep;
        by primary_fbgn;
        run;

    data &sex._&fru._merged;
        merge &sex._&fru._ind (in=in1) &sex._&fru._rep (in=in2);
        if in1 and in2;
        by primary_fbgn;
        run;

    data &sex._&fru._merged;
        set &sex._&fru._merged;
        rename primary_fbgn=FBgn_cat;
        run;

    proc sort data=&sex._&fru._merged;
        by fbgn_cat;
        run;

%mend merge_ind_rep;

%merge_ind_rep(female,fruA);
%merge_ind_rep(female,fruB);
%merge_ind_rep(female,fruC);

%merge_ind_rep(male,fruA);
%merge_ind_rep(male,fruB);
%merge_ind_rep(male,fruC);

/* */
proc sort data=FRU.results;
    by fbgn_cat;
    run;

%macro merge_results(sex,fru);
data &sex._&fru._indrep_results;
    merge &sex._&fru._merged (in=in1) FRU.results (in=in2);
    by fbgn_cat;
    if in1 and in2;
    run;
%mend merge_results;
%merge_results(female,frua);
%merge_results(female,frub);
%merge_results(female,fruc);
%merge_results(male,frua);
%merge_results(male,frub);
%merge_results(male,fruc);

/* Look at Constitutive Common Alternative exons. Because of the naming scheme,
 * it was easier to copy and paste then to use a macro. I used Michelle's
 * significance scheme (FDR flag for a given contrast was significant) and the
 * difference in mean rpkm was greater then 2-fold.
 */

/* Female Fru A contrast 17 
 * 2 FBGN's that were Induced and repressed
 * FBGN0052582
 *      All induced exons were common
 *      All repressed exons were constitutive
 * FBGN0086613
 *      All induced exons were constitutive
 *      All repressed exons were constitutive
 */
data Female_frua_indrep_results2;
    set Female_frua_indrep_results;
    fold = AH_CSFemale_mean_rpkm / AH_Female_FruM_A__mean_rpkm;
    if fold eq '' then delete;
    if fold ge 2 then flag_rep = 1;
    else flag_rep = 0;
    if fold le .5 then flag_ind = 1;
    else flag_ind = 0;
    if fold ge 2 or fold le .5 or flag_fdr_p_contrast_17_20 eq 1 then flag_arb_sig = 1;
    else flag_arb_sig = 0;
    keep fbgn_cat constitutive common alternative flag_rep flag_ind;
    run;

/* Female Fru B contrast 19
 * NO FBGN's were Induced and repressed
 */
data Female_frub_indrep_results2;
    set Female_frub_indrep_results;
    fold = AH_CSFemale_mean_rpkm / AH_Female_FruM_B__mean_rpkm;
    if fold eq '' then delete;
    if fold ge 2 then flag_rep = 1;
    else flag_rep = 0;
    if fold le .5 then flag_ind = 1;
    else flag_ind = 0;
    if fold ge 2 or fold le .5 or flag_fdr_p_contrast_19_20 eq 1 then flag_arb_sig = 1;
    else flag_arb_sig = 0;
    keep fbgn_cat constitutive common alternative flag_rep flag_ind;
    run;

/* Female Fru C contrast 21 
 * No FBGN's were Induced and Repressed
 */
data Female_fruc_indrep_results2;
    set Female_fruc_indrep_results;
    fold = AH_CSFemale_mean_rpkm / AH_Female_FruM_C__mean_rpkm;
    if fold eq '' then delete;
    if fold ge 2 then flag_rep = 1;
    else flag_rep = 0;
    if fold le .5 then flag_ind = 1;
    else flag_ind = 0;
    if fold ge 2 or fold le .5 or flag_fdr_p_contrast_21_20 eq 1 then flag_arb_sig = 1;
    else flag_arb_sig = 0;
    keep fbgn_cat constitutive common alternative flag_rep flag_ind;
    run;

/* Male Fru A contrast 11 
 * 3 FBGN's were induced and repressed
 * FBGN0031993
 *      Induced exons were constitutive or alternative
 *      All repressed exons were alternative
 * FBGN0033159
 *      Induced exons were constitutive or alternative
 *      All repressed exons were alternative
 * FBGN0037970
 *      All Induced exons were constitutive 
 *      All repressed exons were constitutive
 */
data male_frua_indrep_results2;
    set male_frua_indrep_results;
    fold = AH_CS_mean_rpkm / AH_Male_FruM_A__mean_rpkm;
    if fold eq '' then delete;
    if fold ge 2 then flag_rep = 1;
    else flag_rep = 0;
    if fold le .5 then flag_ind = 1;
    else flag_ind = 0;
    if fold ge 2 or fold le .5 or flag_fdr_p_contrast_11_20 eq 1 then flag_arb_sig = 1;
    else flag_arb_sig = 0;
    keep fbgn_cat constitutive common alternative flag_rep flag_ind;
    run;

/* Female Fru B contrast 13
 * No FBGN's were Induced and Repressed
 */
data male_frub_indrep_results2;
    set male_frub_indrep_results;
    fold = AH_CS_mean_rpkm / AH_Male_FruM_B__mean_rpkm;
    if fold eq '' then delete;
    if fold ge 2 then flag_rep = 1;
    else flag_rep = 0;
    if fold le .5 then flag_ind = 1;
    else flag_ind = 0;
    if fold ge 2 or fold le .5 or flag_fdr_p_contrast_13_20 eq 1 then flag_arb_sig = 1;
    else flag_arb_sig = 0;
    keep fbgn_cat constitutive common alternative flag_rep flag_ind;
    run;

/* Female Fru C contrast 15 
 * 7 FBGN's were induced and repressed
 * FBGN0031016
 *      Induced exons were common or constitutive
 *      All repressed exons were constitutive
 * FBGN0032470
 *      All Induced exons were constitutive
 *      All repressed exons were constitutive
 * FBGN0039164
 *      All Induced exons were constitutive
 *      All repressed exons were constitutive
 * FBGN0039632
 *      All Induced exons were constitutive
 *      All repressed exons were constitutive
 * FBGN0052506
 *      All Induced exons were constitutive
 *      All repressed exons were constitutive
 * FBGN0259213
 *      All Induced exons were constitutive
 *      All repressed exons were constitutive
 * FBGN0261574
 *      All Induced exons were constitutive
 *      All repressed exons were constitutive
 */
data male_fruc_indrep_results2;
    set male_fruc_indrep_results;
    fold = AH_CS_mean_rpkm / AH_Male_FruM_C__mean_rpkm;
    if fold eq '' then delete;
    if fold ge 2 then flag_rep = 1;
    else flag_rep = 0;
    if fold le .5 then flag_ind = 1;
    else flag_ind = 0;
    if fold ge 2 or fold le .5 or flag_fdr_p_contrast_15_20 eq 1 then flag_arb_sig = 1;
    else flag_arb_sig = 0;
    keep fbgn_cat constitutive common alternative flag_rep flag_ind;
    run;

