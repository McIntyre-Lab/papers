/*******************************************************************************
* Filename: misc_test_cegsV_inr_corr_with_fusion.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Look at if there is a single fusion that is driving the
* relationship between InR and Sxl. See how the dsx relationships change.
*
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname cegs '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

data mydat;
    set SEM.inr_gene_and_fusion;
    drop line;
    run;

ods listing close;
ods html  path='!MCLAB/cegs_sem_sd_paper/analysis_output/splicing' file='correlations.html';
    proc corr data=mydat;
        var : ;
        run;
ods html close;
ods listing;

