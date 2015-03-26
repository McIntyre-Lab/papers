/********************************************************************************
* Michelle has been giving gene lists for genes that are differentially
* expressed according to her critera. For the Fru Null datasets, we have some
* differences. I am going to apply her filter and see exactly what is the cause
* of this.
*
* Criteria:
*       FDR >0.2 for all 4 fru null comparisons
*       Average Fold difference > 1 for each of the fru null comparisons
********************************************************************************/

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';

data mydata;
    set fru.results_plus_gov2;
    run;

/* Create A subset that removes un-analyzed fusions and multigene fusions */
    data subset;
        set mydata;
        if flag_fdr_p_contrast_7_20 = ' ' then delete;
        if flag_fdr_p_contrast_8_20  = ' ' then delete;
        if flag_fdr_p_contrast_9_20  = ' ' then delete;
        if flag_fdr_p_contrast_10_20 = ' ' then delete;
        where FBgn_cat ^? ';' ;
        *keep fusion_id FBgn_cat AH_BerM_mean_rpkm AH_CS_mean_rpkm AH_FruP14_440_mean_rpkm AH_FruW12_ChaM5_mean_rpkm flag_fdr_p_contrast_7_20 flag_fdr_p_contrast_8_20 flag_fdr_p_contrast_9_20 flag_fdr_p_contrast_10_20;
        run;

/* Identify Fusions that have sig FDR for all comparisons */
    data flag_fdr_four;
        set subset;
        if flag_fdr_p_contrast_7_20 = 1 and flag_fdr_p_contrast_8_20 = 1 and flag_fdr_p_contrast_9_20 = 1 and flag_fdr_p_contrast_10_20 = 1 then flag_fdr_four = 1; else flag_fdr_four = 0;
        if mean(fdr_p_contrast_7,fdr_p_contrast_8,fdr_p_contrast_9,fdr_p_contrast_10) < .2 then flag_fdr_mean_four = 1; else flag_fdr_mean_four = 0;
        keep fusion_id flag_fdr_four flag_fdr_mean_four;
        run;



/* Identify Fusions that have sig FOLD CHANGE for all comparisons */
    data fold;
        set subset;
        fold_P14CS =  AH_CS_mean_rpkm / AH_FruP14_440_mean_rpkm;
        fold_P14BER = AH_BerM_mean_rpkm / AH_FruP14_440_mean_rpkm;
        fold_C5CS =   AH_CS_mean_rpkm / AH_FruW12_ChaM5_mean_rpkm;
        fold_C5BER =  AH_BerM_mean_rpkm / AH_FruW12_ChaM5_mean_rpkm;
        keep Fusion_id fold_P14CS fold_P14BER fold_C5CS fold_C5BER;
        run;

    data flag_up_down;
        set fold;
        if mean(fold_P14CS,fold_P14BER,fold_C5CS,fold_C5BER) gt 1 then flag_up = 1; else flag_up = 0;
        if mean(fold_P14CS,fold_P14BER,fold_C5CS,fold_C5BER) lt 1 then flag_down = 1; else flag_down = 0;
        keep fusion_id flag_up flag_down;
        run;

/* Make a Combined dataset */
proc sort data=flag_fdr_four;
    by fusion_id;
    run;

proc sort data=flag_up_down;
    by fusion_id;
    run;

data merged1;
    merge flag_fdr_four flag_up_down;
    by fusion_id;
    run;

/* Create induced list */
    data fusion_induced_repressed;
        set merged1;
        if flag_fdr_four = 1 and flag_up = 1 then flag_induced = 1; else flag_induced = 0;
        if flag_fdr_mean_four = 1 and flag_up = 1  then flag_induced2 = 1; else flag_induced2 = 0;
        if flag_fdr_four = 1 and flag_down = 1 then flag_repressed = 1; else flag_repressed = 0;
        if flag_fdr_mean_four = 1 and flag_down = 1 then flag_repressed2 = 1; else flag_repressed2 = 0;
        keep fusion_id flag_induced flag_induced2 flag_repressed flag_repressed2;
        run;

    proc freq data = fusion_induced_repressed;
        table flag_induced;
        table flag_induced2;
        table flag_repressed;
        table flag_repressed2;
        run;

/* Merge on FBgn_cat and collapse to gene level */
proc sort data=subset;
    by fusion_id;
    run;

proc sort data=fusion_induced_repressed;
    by fusion_id;
    run;

data merged2;
    merge subset (in=in1) fusion_induced_repressed (in=in2);
    by fusion_id;
    if in1 and not in2 then flag_induced = 0 and flag_induced2 = 0 and flag_repressed = 0 and flag_repressed2 = 0;
    keep fusion_id FBgn_cat flag_induced flag_induced2 flag_repressed flag_repressed2;
    run;

proc sort data=merged2;
    by FBgn_cat;
    run;

proc means data=merged2 noprint;
    by FBgn_cat;
    output out=sum sum(flag_induced)= sum(flag_induced2)= sum(flag_repressed)= sum(flag_repressed2)= /autoname;
    run;

proc freq data=sum;
    table flag_induced_sum;
    table flag_induced2_sum;
    table flag_repressed_sum;
    table flag_repressed2_sum;
    run;

data null_induced;
    set sum;
    where flag_induced_sum > 0;
    keep FBgn_cat;
    run;

data null_induced2;
    set sum;
    where flag_induced2_sum > 0;
    keep FBgn_cat;
    run;

data null_repressed;
    set sum;
    where flag_repressed_sum > 0;
    keep FBgn_cat;
    run;

data null_repressed2;
    set sum;
    where flag_repressed2_sum > 0;
    keep FBgn_cat;
    run;

proc export data=null_induced
            outfile= '/home/jfear/mclab/arbeitman/arbeitman_fru_network/exported_data_from_michelle/Induced_Fru_m_null_jmf.tab'
            dbms=CSV replace;
            putnames=no;
            run;

proc export data=null_repressed
            outfile= '/home/jfear/mclab/arbeitman/arbeitman_fru_network/exported_data_from_michelle/Repressed_Fru_m_null_jmf.tab'
            dbms=CSV replace;
            putnames=no;
            run;

data FRU.Null_male_fusions_ind_rep;
    set merged2;
    rename flag_induced = flag_null_male_ind;
    rename flag_repressed = flag_null_male_rep;
    keep fusion_id flag_induced flag_repressed;
    run;
