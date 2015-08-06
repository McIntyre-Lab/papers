/*******************************************************************************
* Filename: splicing_InR_splicing_model.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Run a set of linear models to decide which fusions and genotypes
* to look at.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname cegs '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

proc sort data=sem.cegsV_splice_data;
    by symbol_cat fusion_id;
    run;

proc glm data=sem.cegsV_splice_data ;
    by symbol_cat fusion_id;
    class line;
    model uq_log_uq_center=line ;
    ods output  modelanova=glmbyfusion;
    run;
    quit;
