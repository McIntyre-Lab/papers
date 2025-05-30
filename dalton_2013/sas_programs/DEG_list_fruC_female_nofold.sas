/********************************************************************************
* Michelle has been giving gene lists for genes that are differentially
* expressed according to her critera. For the Fru Null datasets, we have some
* differences. I am going to apply her filter and see exactly what is the cause
* of this.
*
* Criteria:
*       FDR <0.2 for all 4 fru null comparisons
*       Directions the same for each of the fru comparisons
********************************************************************************/

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';

data mydata;
    set fru.results_plus_gov2;
    run;

/* Create A subset that removes un-analyzed fusions and multigene fusions */
    data subset;
        set mydata;
        if flag_fdr_p_contrast_21_20 = ' ' then delete;
        if flag_fdr_p_contrast_22_20  = ' ' then delete;
        where FBgn_cat ^? ';';
        run;

/* Identify Fusions that have sig FDR for all comparisons */
    data flag_fdr_all;
        set subset;
        if flag_fdr_p_contrast_21_20 = 1 and flag_fdr_p_contrast_22_20 = 1 then flag_fdr_all = 1; else flag_fdr_all= 0;
        keep fusion_id flag_fdr_all;
        run;


/* Identify Fusions direction of change for all comparisons */
    data fold;
        set subset;
        fold_CS =  AH_CSFemale_mean_rpkm / AH_FEMALE_FruM_C__mean_rpkm;
        fold_BER = AH_BerF_mean_rpkm / AH_FEMALE_FruM_C__mean_rpkm;
        keep Fusion_id fold_CS fold_BER;
        run;

    data flag_up_down;
        set fold;
        if fold_CS lt 1 and fold_BER lt 1 then flag_up = 1; else flag_up = 0;
        if fold_CS gt 1 and fold_BER gt 1 then flag_down = 1; else flag_down = 0;
        keep fusion_id flag_up flag_down;
        run;

/* Make a Combined dataset */
proc sort data=flag_fdr_all;
    by fusion_id;
    run;

proc sort data=flag_up_down;
    by fusion_id;
    run;

data merged1;
    merge flag_fdr_all flag_up_down;
    by fusion_id;
    run;

/* Create induced list */
    data fusion_induced_repressed;
        set merged1;
        if flag_fdr_all = 1 and flag_up = 1 then flag_induced = 1; else flag_induced = 0;
        if flag_fdr_all = 1 and flag_down = 1 then flag_repressed = 1; else flag_repressed = 0;
        keep fusion_id flag_induced flag_repressed ;
        run;

    proc freq data = fusion_induced_repressed;
        table flag_induced;
        table flag_repressed;
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
    if in1 and not in2 then flag_induced = 0 and flag_repressed = 0;
    keep fusion_id FBgn_cat flag_induced flag_repressed ;
    run;

proc sort data=merged2;
    by FBgn_cat;
    run;

proc means data=merged2 noprint;
    by FBgn_cat;
    output out=sum sum(flag_induced)= sum(flag_repressed)= /autoname;
    run;

proc freq data=sum;
    table flag_induced_sum;
    table flag_repressed_sum;
    run;

data induced;
    set sum;
    where flag_induced_sum > 0;
    keep FBgn_cat;
    run;

data repressed;
    set sum;
    where flag_repressed_sum > 0;
    keep FBgn_cat;
    run;

proc export data=induced
            outfile= '/home/jfear/mclab/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_female_ind_1fold_jmf.tab'
            dbms=CSV replace;
            putnames=no;
            run;

proc export data=repressed
            outfile= '/home/jfear/mclab/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_female_rep_1fold_jmf.tab'
            dbms=CSV replace;
            putnames=no;
            run;

data FRU.fruC_female_fus_ind_rep_nfold;
    set merged2;
    rename flag_induced = flag_fruC_female_ind_nofold;
    rename flag_repressed = flag_fruC_female_rep_nofold;
    keep fusion_id flag_induced flag_repressed;
    run;

proc export data=FRU.fruC_female_fus_ind_rep_nfold
            outfile= '/home/jfear/mclab/arbeitman/arbeitman_fru_network/reports_external/FruM_C_female_flag_ind_rep_nofold.csv'
            dbms=CSV replace;
            putnames=yes;
            run;
