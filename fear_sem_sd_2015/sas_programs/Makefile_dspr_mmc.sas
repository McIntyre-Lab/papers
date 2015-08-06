libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Modulated Modularity Clustering on Collapsed Isoforms */
    * Prepare datset for upload to the MMC website. Genes should be rows and
    * observations should be columns.
    *
    * INPUT: SEM.dsrp_sbs_combine_sym
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/mmc/dsrp_isoform_sex_det_for_mmc.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_export_for_MMC.sas';

/* Modulated Modularity Clustering on All Isoforms */
    * Prepare datset for upload to the MMC website. Genes should be rows and
    * observations should be columns.
    *
    * INPUT: SEM.dsrp_sbs_symbol
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/mmc/dsrp_all_isoform_sex_det_for_mmc.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_export_for_MMC_all_iso.sas';

/* Modulated Modulated on Sex Det Genes */
    data tmp;
        set SEM.dsrp_sbs_gene_level_sym;
        line = strip(patRIL) || '_' || strip(matRIL);
        drop CG: patRIL matRIL;
        run;

    proc transpose data=tmp out=flip;
        id line;
        run;

        
    proc export data=flip outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/mmc/dsrp_sex_det_genes_for_mmc.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Clean up;
    proc datasets nolist;
        delete tmp flip;
        run; quit;
