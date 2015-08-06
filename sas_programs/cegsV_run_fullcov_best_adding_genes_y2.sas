/*******************************************************************************
* Filename: cegsV_run_fullcov_best_adding_genes_y2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

proc sort data=SEM.cegsV_ag_yp2_stack_bic;
    by bic;
    run; *Best model FBgn0028325 Model 33, BIC = 436.04;

proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
lismod
yvar = Sxl Yp2 dsx fru tra FBgn0028325, 
xvar = Spf45 fl_2_d her ix snf tra2 vir
;
matrix _BETA_ 
[1,6] = beta1,
[2,3] = beta2,
[3,5] = beta3,
[4,5] = beta4,
[5,1] = beta5
;
matrix _GAMMA_ 
[1,1] = gamma1,
[1,2] = gamma2,
[1,7] = gamma3,
[2,3] = gamma4,
[2,4] = gamma5,
[3,6] = gamma6,
[4,6] = gamma7,
[5,2] = gamma8,
[5,7] = gamma9,
[6,5] = gamma10
;
run; *BIC = 362.8397;
