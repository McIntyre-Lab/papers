libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Gaussian Graphic model on ALL ISOFORMS in DSRP */
    * Create and export dataset for gaussian graphical models using the sex
    * determination subset. Then run GGM R script. GGM allows you to select
    * edges by an FDR cutoff or by specifiying a number of edges. Two
    * different FDR cutoffs were used (0.2, 0.8) along with outputing the
    * top 20 edges. GGM was run on the entire DSRP dataset.
    * 
    * INPUT: SEM.dsrp_stack
    *
    * RSCRIPT: !MCLAB/cegs_sem_sd_paper/r_programs/dsrp_ggm_all_isoforms.R
    *
    * FILES: $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.dot
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.png 
    *
    * NOTE:  I could not get sas to export all (11065) columns. I had to
    * export the intermediate csv file using JMP genomics.
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_create_ggm_dataset_DSRP.sas';

/* Gaussian Graphic model on GENES in DSRP */
    * Create and export dataset for gaussian graphical models using the sex
    * determination subset. Then run GGM R script. GGM allows you to select
    * edges by an FDR cutoff or by specifiying a number of edges.  An FDR
    * cutoff of 0.2 was used. GGM was run on the entire DSRP dataset at the
    * gene level.
    * 
    * INPUT: SEM.dsrp_sbs_gene_level_sym 
    *
    * RSCRIPT: !MCLAB/cegs_sem_sd_paper/r_programs/dsrp_ggm_gene.R
    *
    * FILES: $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_gene_FDR2.dot
    *        /home/jfear/tmp/dsrp_by_gene_sbs.txt
    *
    * NOTE:  I could not get sas to export all (7424) columns. I had to
    * export the intermediate csv file using JMP genomics.
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_create_gene_level_ggm_dataset.sas';

/* Gaussian Graphic Model on sex det subset */
    * Create and export dataset for gaussian graphical models using the sex
    * determination subset. Then run GGM R script. GGM allows you to select
    * edges by an FDR cutoff or by specifiying a number of edges. Two
    * different FDR cutoffs were used (0.2, 0.8) along with outputing the
    * top 20 edges. Also Isoforms and gene level variables are used.
    *
    * INPUT: SEM.dsrp_sex_det_sbs_combine_sym;
    *        SEM.dsrp_sex_det_sbs_gene_level_sym;
    * 
    * RSCRIPTS: !MCLAB/cegs_sem_sd_paper/r_programs/dsrp_ggm_sex_det_subset_gene_v2.R
    *           !MCLAB/cegs_sem_sd_paper/r_programs/dsrp_ggm_sex_det_subset_v2.R
    *
    * FILES: $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_isoforms_FDR2.dot
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_isoforms_FDR2.png 
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_isoforms_TOP20.dot
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_isoforms_TOP20.png 
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_gene_FDR2.dot
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_gene_FDR2.png 
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_gene_TOP20.dot
    *        $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_sex_det_gene_TOP20.png 
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_create_ggm_dataset_sex_det_subset.sas';

/* Coverage vs Connectedness */
    * Trying to address the question of coverage vs connectedness. I need to
    * plot the number of neighbors in relationship to averge expression level.
    *
    * INPUT: SEM.DSRP_stack
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_neighborhood_analysis_v3.csv
    * 
    * FILE: MCLAB/cegs_sem_sd_paper/analysis_output/ggm/expression_v_neighbohood.png
    *       MCLAB/cegs_sem_sd_paper/analysis_output/ggm/expression_v_std.png
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_connectedness_v_expression.sas';

/* Create GGN gene list */
    * I used ggm_neighborhood_analysis_v3.py to generate a table of counts for
    * genes in the 1-step neighborhood and genes in the 2-step
    * neighborhood. This table also includes which genes in the sex
    * hierarchy are in the primary and secondary neighborhoods. Create a
    * gene list, where a gene needs to have 2 or more genes from the sex
    * hierarchy in thier primary neighborhood.
    *
    * INPUT: !MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_neighborhood_analysis_v3.csv
    *
    * DATASET: SEM.dsrp_ggn_primary_neigh_gene_list;
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_primary_neighborhood_gene_list.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_ggm_neighborhood_enrichment.sas';
