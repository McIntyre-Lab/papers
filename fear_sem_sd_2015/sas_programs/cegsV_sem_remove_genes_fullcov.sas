/*******************************************************************************
* Filename: cegsV_sem_remove_genes.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Remove each gene_droppedin the sex hierarchy one at time and see how
* that affects BIC
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
lismod
yvar = Sxl Yp2 dsx fru tra, 
xvar = Spf45 fl_2_d her ix snf tra2 vir
;
matrix _BETA_ 
[2,3] = beta1,
[3,5] = beta2,
[4,5] = beta3,
[5,1] = beta4
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
run;

/* Baseline */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'baseline';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* Spf45 */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx her ix -> Yp2 
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'Spf45';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* fl_2_d */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        Spf45 snf vir -> Sxl,
        vir Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx her ix -> Yp2 
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'fl_2_d';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* snf */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'snf';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* vir */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf -> Sxl,
        fl_2_d Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'vir';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* Sxl */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> tra,
        tra2 tra -> fru dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'Sxl';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* tra */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        tra2 Sxl -> fru dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'tra';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* tra2 */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra -> fru dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'tra2';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* fru */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> dsx,
        dsx her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'fru';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* dsx */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru,
        her ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'dsx';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* her */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx ix -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'her';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* ix */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru dsx,
        dsx her -> Yp2
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'ix';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* Yp2 */
proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    path
        fl_2_d Spf45 snf vir -> Sxl,
        fl_2_d vir Sxl -> tra,
        tra2 tra -> fru dsx
    ;
    run;

data bic;
    length gene_dropped $10.;
    retain gene_dropped fitvalue;
    set fitstat;
    where indexCode eq 312;
    gene_dropped = 'Yp2';
    keep gene_dropped fitvalue;
    run;

proc append data=bic base=combined;
run;

/* Sort by BIC */
data combined2;
    set combined;
    rename fitvalue = bic;
    run;

proc sort data=combined2;
    by bic;
    run;

/* Save permanent dataset */
data SEM.cegsV_sem_remove_genes_fullcov;
    set combined2;
    run;

/* Clean up */
proc datasets ;
    delete bic;
    delete fitstat;
    delete combined;
    delete combined2;
    run; quit;

