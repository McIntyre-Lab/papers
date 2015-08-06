libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel530 "!MCLAB/useful_dmel_data/flybase530/sas_data";
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);


/* Create Gene List for Male vs Female controls */
    * Using the CS and Berlin control data from Michelle, create a gene list.
    * Require that the contrast to both controls was significant with an FDR
    * 0.05 and that the direction was the same. 
    *
    * INPUT: dmel530.Fb530_si_fusions_unique_flagged
    *        dsx.exp_means
    *        dsx.flag_fdr_contrast_by_fusion
    *
    * DATASET: SEM.dsx_ctrl_female_induced      246 genes
    *          SEM.dsx_ctrl_female_repressed    205 genes
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsx_create_dsx_male_v_female_dataset.sas';

/* Create Gene List for dsxD females */
    * Using the dsxD data from Michelle, create a gene list. Require that the
    * contrast to both controls was significant with an FDR 0.05 and that the
    * direction was the same. 
    *
    * INPUT: dmel530.Fb530_si_fusions_unique_flagged
    *        dsx.exp_means
    *        dsx.flag_fdr_contrast_by_fusion
    *
    * DATASET: SEM.dsxD_induced                 238 genes 
    *          SEM.dsxD_repressed               239 genes
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsx_create_dsxD_dataset.sas';

/* Create Gene List for dsxNULL females */
    * Using the dsxNULL female data from Michelle, create a gene list. Require
    * that the contrast to both controls was significant with an FDR 0.05 and
    * that the direction was the same. 
    *
    * INPUT: dmel530.Fb530_si_fusions_unique_flagged
    *        dsx.exp_means
    *        dsx.flag_fdr_contrast_by_fusion
    *
    * DATASET: SEM.dsxNullF_induced             340 genes
    *          SEM.dsxNullF_repressed           208 genes
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsx_create_dsxNULL_female_dataset.sas';

/* Merge Gene lists together */
    * Now that I have distinct gene lists I can merge them together and create
    * a Venn diagram to see how they overlapp.
    *
    * INPUT: SEM.dsxD_induced                 238 genes 
    *        SEM.dsxD_repressed               239 genes
    *        SEM.dsxNullF_induced             340 genes
    *        SEM.dsxNullF_repressed           208 genes
    *        SEM.dsx_ctrl_female_induced      246 genes
    *        SEM.dsx_ctrl_female_repressed    205 genes
    *        dmel530.Fb530_si_fusions_unique_flagged
    *
    * Rscript: $MCLAB/cegs_sem_sd_paper/r_programs/dsx_venn.R
    *
    * FILES: $MCLAB/cegs_sem_sd_paper/analysis_output/dsx/induced_venn.png
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/dsx/repressed_venn.png
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/dsx/both_venn.png
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/dsx/induced_venn.svg
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/dsx/repressed_venn.svg
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/dsx/both_venn.svg
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsx_merge_gene_lists.sas';

data induced repressed;
    set SEM.validation_set_bic12 ;
    if flag_dsxNullF_induced = 1 then do;
    description = 'dsxNullF Induced';
    output induced;
    end;
    if flag_dsxNullF_repressed = 1 then do;
    description = 'dsxNullF Repressed';
    output repressed;
    end;
    keep primary_fbgn description;
    run;

proc sort data=induced nodupkey;
    by primary_fbgn;
    run;

proc sort data=repressed nodupkey;
    by primary_fbgn;
    run;

proc export data=induced outfile='/home/jfear/mclab/useful_dmel_data/gene_lists/exported_data/dsx_null_induced.csv' dbms=csv replace;
    putnames=yes;
    run;

proc export data=repressed outfile='/home/jfear/mclab/useful_dmel_data/gene_lists/exported_data/dsx_null_repressed.csv' dbms=csv replace;
    putnames=yes;
    run;
