libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Modulated Modularity Clustering */
    * Prepare datset for upload to the MMC website. Genes should be rows and
    * observations should be columns.
    *
    * INPUT: SEM.cegsV_sex_det_stack
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/mmc/cegsV_sex_det_data_for_mmc.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_export_for_MMC.sas';

/* Modulated modularity cluster sex det genes */
    proc transpose data=SEM.cegsV_by_gene_sex_det_sbs  out=flip;
        id line;
        run;
        
    proc export data=flip outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/mmc/cegsV_sex_det_gene_for_mmc.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Clean up;
    proc datasets nolist;
        delete flip;
        run; quit;
