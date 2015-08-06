/*******************************************************************************
* Filename: cegsV_run_fullcov_best_adding_links_yp2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries

*/
proc sort data=SEM.cegsV_al_yp2_stack_bic;
    by bic;
    run; * model 3 Sxl -> Fru BIC = 483.73586;

proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
lismod
yvar = Sxl Yp2 dsx fru tra, 
xvar = Spf45 fl_2_d her ix snf tra2 vir
;
matrix _BETA_ 
[2,3] = beta1,
[3,5] = beta2,
[4,1] = beta3,
[4,5] = beta4,
[5,1] = beta5
;
matrix _GAMMA_ 
[1,1] = gamma1,
[1,2] = gamma2,
[1,5] = gamma3,
[1,7] = gamma4,
[2,3] = gamma5,
[2,4] = gamma6,
[3,6] = gamma7,
[4,6] = gamma8,
[5,2] = gamma9,
[5,7] = gamma10
;
run; * 414.8469 ;
