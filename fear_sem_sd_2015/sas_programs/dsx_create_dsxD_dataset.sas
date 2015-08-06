/********************************************************************************
* I am wanting to utilize the dsxD data to identify potential genes affected
* downstream of DSX. The dsxD is the male isoform that was expressed in the
* female. So genes disregulated by dsxD are downstream of dsx.
********************************************************************************/
/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel530 "!MCLAB/useful_dmel_data/flybase530/sas_data";
    libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Combine means and fdr flags */
    proc sort data=dmel530.Fb530_si_fusions_unique_flagged;
        by fusion_id;
        run;

    proc sort data=dsx.exp_means;
        by fusion_id;
        run;

    proc sort data=dsx.flag_fdr_contrast_by_fusion;
        by fusion_id;
        run;

    data mydata;
        merge dmel.Fb530_si_fusions_unique_flagged (in=in1) dsx.exp_means (in=in2) 
              dsx.flag_fdr_contrast_by_fusion (in=in3);
        by fusion_id;
        run;

/* Create A subset that removes un-analyzed fusions and multigene fusions */
    data subset;
        set mydata;
        if flag_fdr_p_contrast_10_05 = ' ' then delete;
        if flag_fdr_p_contrast_8_05  = ' ' then delete;
        where FBgn_cat ^? ';';
        run;

/* Identify Fusions that have sig FDR for all comparisons */
    * BerF-dsxD(10), CSF-dsxD(8)
    ;
    data flag_fdr_all;
        set subset;
        if flag_fdr_p_contrast_8_05 = 1 and flag_fdr_p_contrast_10_05 = 1 then flag_fdr_all = 1; else flag_fdr_all= 0;
        keep fusion_id flag_fdr_all;
        run;

/* Identify direction comparisons */
    data direction;
        set subset;
        dir_BER =  AH_BerF_mean_rpkm / AH_dsxD_mean_rpkm;
        dir_CS = AH_CSFemale_mean_rpkm / AH_dsxD_mean_rpkm;
        keep Fusion_id dir_CS dir_BER;
        run;

    data flag_up_down;
        set direction;
        if dir_CS lt 1 and dir_BER lt 1 then flag_up = 1; else flag_up = 0;
        if dir_CS gt 1 and dir_BER ge 1 then flag_down = 1; else flag_down = 0;
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

    data SEM.dsxD_induced;
        set sum;
        where flag_induced_sum > 0;
        keep FBgn_cat;
        run;    
        * FDR .20 = 625 obs;
        * FDR .10 = 370 obs;
        * FDR .05 = 238 obs;

    data SEM.dsxD_repressed;
        set sum;
        where flag_repressed_sum > 0;
        keep FBgn_cat;
        run;
        * FDR .20 = 834 obs;
        * FDR .10 = 433 obs;
        * FDR .05 = 239 obs;

/* Clean-up */
proc datasets;
    delete direction;
    delete flag_fdr_all;
    delete flag_up_down;
    delete fusion_induced_repressed;
    delete merged1;
    delete merged2;
    delete mydata;
    delete subset;
    delete sum;
    run; quit;
