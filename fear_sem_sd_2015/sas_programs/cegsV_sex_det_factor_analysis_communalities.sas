/********************************************************************************
* Look at communalities by running several factor analysis iteratively.
********************************************************************************/
/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Run Factor Analysis All Fusions */
    data sex_det;
        set SEM.cegsV_by_fusion_sex_det_sbs;
        drop line;
        run;

    proc transpose data=sex_det(obs=1) out=flip;
        var :;
        run;

    data cegsV_factor_communal;
        set flip;
        drop col1;
        run;
        
    proc sort data=cegsV_factor_communal;
        by _name_;
        run;
    
    %macro iter_nfact(num);
        proc printto LOG="!MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_communal/cegsV_sex_det_factor_analysis_&num..log" PRINT="!MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_communal/cegsV_sex_det_factor_analysis_&num..lst" NEW;
            run;
        proc factor data=sex_det simple method=PRINCIPAL scree rotate=varimax nfact=&num round flag=.4 outstat=factor;
            var :;
            run;
        proc printto LOG=LOG PRINT=PRINT;
        run;

        data com;
            set factor;
            if _type_ eq 'COMMUNAL';
            drop _type_ _name_;
            run;

        proc transpose data=com out=flip;
            var:;
            run;

        proc sort data=flip;
            by _name_;
            run;

        data merged;
            merge cegsV_factor_communal (in=in1) flip (in=in2);
            by _name_;
            rename col1 = num_factor_&num;
            run;

        data cegsV_factor_communal;
            set merged;
            run;
    %mend;
    %iter_nfact(2);
    %iter_nfact(3);
    %iter_nfact(4);
    %iter_nfact(5);
    %iter_nfact(6);
    %iter_nfact(7);
    %iter_nfact(8);
    %iter_nfact(9);
    %iter_nfact(10);
    %iter_nfact(11);
    %iter_nfact(12);
    %iter_nfact(13);
    %iter_nfact(14);
    %iter_nfact(15);
    %iter_nfact(16);
    %iter_nfact(17);
    %iter_nfact(18);
    %iter_nfact(19);
    %iter_nfact(20);
    %iter_nfact(21);
    %iter_nfact(22);
    %iter_nfact(23);
    %iter_nfact(24);
    %iter_nfact(25);
    %iter_nfact(26);
    %iter_nfact(27);
    %iter_nfact(28);
    %iter_nfact(29);
    %iter_nfact(30);

    data SEM.cegsV_sex_det_factor_commun;
        set cegsV_factor_communal;
        run;

/* Clean up */
    proc datasets nolist;
        delete com;
        delete cegsV_factor_communal;
        delete factor;
        delete flip;
        delete merged;
        delete sex_det;
        run; quit;
