/*******************************************************************************
* Filename: splicing_calc_inr_ratios.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Make table of fusion ratios. I will use this data to make cell
* plots. I am specifically interested in the ratio of S57607/S57601 and
* S57607/S57609. 
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname cegs '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data SEM.inr_fusion_ratio;
    set SEM.inr_gene_and_fusion;
    S57601_SI_S57607_SI = S57601_SI/S57607_SI ;
    F57609_SI_S57607_SI = F57609_SI/S57607_SI ;
    keep line S57601_SI_S57607_SI F57609_SI_S57607_SI S57601_SI S57607_SI F57609_SI;
    run;
