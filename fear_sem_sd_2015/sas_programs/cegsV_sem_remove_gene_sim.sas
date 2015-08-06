/*******************************************************************************
* Filename: cegsV_sem_remove_genes.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Remove each gene_droppedin the sex hierarchy one at time and see how
* that affects BIC
*
*******************************************************************************/

/* Split sysparm */
    data _null_;
        length sysparm express param value $200;
        sysparm = symget('sysparm');
        do i=1 to 2;
            express = left(scan(sysparm, i, ','));
            param = left(scan(express, 1, '='));
            value = left(scan(express, 2, '='));
            call symput(param, trim(left(value)));
        end;
        run;

/* Libraries */
    libname sem "&lib";

proc import datafile="&mydat" out=SEM.mySimDat dbms=csv replace;
    getnames=yes;
    run;

/* Baseline */
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> vir =0,
            snf <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            Spf45 snf vir -> Sxl,
            vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf -> Sxl,
            fl_2_d Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> Spf45 =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,

            snf <-> Spf45 =0,
            snf <-> tra2 =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            tra2 Sxl -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,

            Spf45 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru,
            her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx ix -> Yp2 ,

            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    proc calis data=SEM.mySimDat method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,
            snf <-> tra2 =0,

            Spf45 <-> vir =0,
            Spf45 <-> tra2 =0,

            tra2 <-> vir =0
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
    data SEM.fitStat;
        set combined;
        rename fitvalue = bic;
        run;

    proc sort data=SEM.fitStat;
        by bic;
        run;
