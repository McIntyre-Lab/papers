libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Update DSX FBgn Number to FB5.51 */
    * DSX data was aligned to FB5.30, so the FBgns could potentially have
    * changed. Need to check this and update the FBgns so I can merge this to
    * the CEGS adding gene list
    *
    * NOTE: There were 16 primary FBgns in the DSXNull dataset that change in
    * the Fb5.51. One of these FBgn0261701 was not in the FB5.51 secondary_fbgn
    * list, so I had to use cgnumber "CG42737" to relate it to FBgn0262838.
    *
    * Also there are two genes FBgn0005563 and FBgn0040856 that were merged in
    * FB5.51 to FBgn0263111. I only kept FBgn0005563.
    *
    * INPUT: DMEL551.symbol2fbgn
    *        SEM.DSXNullF_repressed
    *        SEM.DSXNullF_induced
    *        !MCLAB/useful_dmel_data/flybase551/flybase_files/fbgn_annotation_ID_fb_2012_06.tsv
    *
    * DATASET: SEM.DSXNullF_repressed_w_fb551
    *          SEM.DSXNullF_induced_w_fb551
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_vs_dsxNull_update_dsx_fbgn.sas';

/* Merge and Flag DSXNull and CEGS Adding Gene List */
    * Compare the gene list from the CEGSV adding genes algorithm to that of
    * the DSXNull list. The cegsV gene list has symbols for the sex det genes
    * and FBgn numbers for everytyhing else. First convert the sex det symbols
    * to FBgn. Then merge the two datasets and make flags.
    *
    * INPUT: SEM.cegsV_ag_nofilter_yp2_gene_list
    *        DMEL551.symbol2fbgn
    *        SEM.Dsxnullf_repressed_w_fb551
    *
    * DATASET: SEM.flag_ag_dsxNull
    * 
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/dsxNullF_repressed_fb551.csv
    *
    * FLAGS: flag_dsxnull_repessed = 1 if in the DSXNullF repressed list
    *        flag_dsxnull_induced = 1 if in the DSXNullF induced list
    *        flag_add = 1 if in the cegsV yp2 adding genes list
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_vs_dsxNull_merge_and_flag.sas';

/* Make gene List and export */
    * Pull out the list of 28 genes that are overlapping in the datast. Merge
    * on gene symbols and output the data.
    * 
    * INPUT: SEM.flag_ag_dsxNull
    *        DMEL551.symbol2fbgn
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_overlap_w_dsxNull_repressed.csv
    *       !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_overlap_w_dsxNull_induced.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_vs_dsxNull_export_gene_list.sas';

/* Compare added genes list to Luo 2011 */
    * Luo 2011 identified 58 genes with dsx binding sites in them. I want
    * to compare their list to the genes added to the GRN.
    *
    * INPUT: SEM.flag_ag_dsxNull
    *        SEM.Luo2011
    *        SEM.cegsV_ag_yp2_flag_ds_dsx
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_w_flag_dsxNull_flag_luo.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_yp2_vs_dsxNull_compare_luo2011.sas';




    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_vs_dsxNull_where_is_INR.sas';
