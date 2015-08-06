libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname norm '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Pull in CEGS and AVG across reps */
    * Grab the virgins from the normalized centered data and average across
    * reps and merge on gene info.
    * 
    * INPUT: CEGS.ccfus_norm_centered
    *
    * DATASET: SEM.cegs_virgin_norm_cent
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_avg_reps.sas';

/* Create Sex Determination Subset */
    * Create a subset of only the sex determination genes.
    *
    * Note: I kept Sex det fusions, even if they were multigene.
    *
    * INPUT: SEM.cegs_virgin_norm_cent
    *
    * DATASET: SEM.cegsV_sex_det_stack
    *
    * Genes Present In Sex Det Dataset: B52 Psi Rbp1 Spf45 Sxl Yp1 Yp2 Yp3
    *                                   dsx fl_2_d fru her ix mub ps snf
    *                                   sqd tra tra2 vir
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_create_sex_det_dataset.sas';

/* QC of Sex Det subset */
    * The input set of fusions should already be pretty good, because they
    * made it through the CEGS normalization process. This is where I
    * removed line and fusions with low coverage. I need to combine fusions
    * to the gene level, but I want to be very stringent with my selection
    * criteria. This is an exploratory script to look at how the fusions in
    * the sex det pathway are behaving. 
    * 
    * INPUT: SEM.cegsV_sex_det_stack
    * 
    * Number of fusions in each sex det gene: 
    *   B52: 7 Psi: 6 Rbp1: 2 Spf45: 3 Sxl: 9 Yp1: 2 Yp2: 2 Yp3: 3 dsx: 4
    *   fl_2_d: 5 fru: 11 her: 3 ix: 1 mub: 12 ps: 8 snf: 2 sqd: 9 tra: 2
    *   tra2: 4 vir: 5
    * 
    * Number of fusions present after droping those that fail normality
    *   B52: 1 Psi: 2 Rbp1: 1 Spf45: 2 Sxl: 1 Yp1: 2 Yp2: 2 Yp3: 3 dsx: 3
    *   fl_2_d: 0 fru: 3 her: 1 ix: 1 mub: 2 ps: 5 snf: 0 sqd: 1 tra: 1
    *   tra2: 0 vir: 2
    * 
    * Looking closer at tra2: 
    *   F40928_SI looks ok inspite of failing the normality cutoff
    * 
    * Looking closer at snf: 
    *   S13633_SI  looks ok inspite of failing the normality cutoff
    *
    * Looking closer at fl_2_d: 
    *   F40441_SI looks ok inspite of failing the normality cutoff
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_qc.sas';

/* Build Gene SBS dataset */
    * Build a gene sbs dataset for use in SEMs. Lauren and I talked and she
    * felt my golden dataset was too stringent, so I am creating a 'bronze'
    * dataset that does not remove fusions that fail normality.
    *
    * (1) Merge on sex det formated symbols
    * (3) Make side-by-side dataset
    * (4) Make a list of genes for iterating over later
    *
    * INPUT: SEM.cegs_virgin_norm_cent
    *        SEM.cegsV_sex_det_stack
    *
    * DATASET: SEM.cegsV_by_gene_stack
    *          SEM.cegsV_by_gene_sbs
    *          SEM.cegsV_by_gene_sex_det_sbs
    *          SEM.cegsV_gene_list
    *          SEM.cegsV_gene2fbgn
    *          SEM.cegsV_by_gene_sbs_mean_center
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/design_file/cegsV_by_gene_sbs_gene_list.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_build_gene_sbs.sas';

/* Build Fusion SBS dataset */
    * Build a fusion sbs dataset.
    *
    * (1) Merge on sex det formated symbols
    * (3) Make side-by-side dataset
    *
    * INPUT: SEM.cegs_virgin_norm_cent
    *        SEM.cegsV_sex_det_stack
    *
    * DATASET: SEM.cegsV_by_fusion_sbs
    *          SEM.cegsV_by_fusion_sex_det_sbs
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_build_fusion_sbs.sas';

/* Build Fusion CNT SBS dataset */
    * Build a fusion sbs dataset, but instead of fusions ID use a counter.
    *
    * (1) Merge on sex det formated symbols
    * (3) Make side-by-side dataset
    *
    * INPUT: SEM.cegs_virgin_norm_cent
    *        SEM.cegsV_sex_det_stack
    *
    * DATASET: SEM.cegsV_by_fusion_cnt_sbs
    *          SEM.cegsV_by_fusion_sex_det_cnt_sbs
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_build_fusion_cnt_sbs.sas';

/* Calculate correlation between isoforms in all genes and combine */
        * Multicollinearity is going to be a large problem. It will be best to
        * collapse isoforms if they have a large correlation coefficient (>=0.75).
        * 
        * Uses a python script to calculate the correlation between isoforms
        * and combine them (AVERAGE) if they have a correlation >=0.75. 
        *
        * Went from 32722 fusions to 30955, not a big change.
        *
        * INPUT: SEM.cegsV_by_fusion_cnt_sbs
        *
        * SCRIPTS: $MCLAB/cegs_sem_sd_paper/scripts/cegs_clusterByCorr.py 
        *
        * DATASET: SEM.cegsV_by_fusion_comb_sbs
        *          
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/cegsV/&isoform..csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/cegsV/&isoform._iso_table.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_calculate_correlation_and_combine.sas';

/* Create Line list */
    * Create a list of lines that we are using in the cegs_sd paper
    * 
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * DATASET: SEM.cegsV_line_list
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/design_file/cegsV_line_list.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_line_list.sas';

/* Randomly Select 1/2 lines */
    * I want to run the CEGS data using half of the genotypes. If our hypothesis
    * about the DSPR is true then we should not add any genes when using half
    * the number of genotypes.
    *
    * INPUT: SEM.cegsV_line_list
    *
    * DATASET: SEM.cegsV_random_subset
    *          SEM.cegsV_by_gene_sbs_subset
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_random_subset.sas';
