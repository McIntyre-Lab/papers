/********************************************************************************
* Using the DSX datasets, I have created a gene list with genes repressed
* when dsxF is removed (dsxNull females).
********************************************************************************/
/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/


/* Run Factor Analysis All Genes */
    proc transpose data=SEM.dsrp_dsx_sbs_symbol(obs=1) out=flip;
        var CG:;
        run;
        
    data factor_all_genes_communal;
        set flip;
        drop col1;
        run;
        
    proc sort data=factor_all_genes_communal;
        by _name_;
        run;
    
    %macro iter_nfact(num);
        proc printto LOG="!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_all_&num..log" PRINT="!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_all_&num..lst" NEW;
            run;
        proc factor data=SEM.dsrp_dsx_sbs_symbol simple method=PRINCIPAL scree rotate=varimax nfact=&num round flag=.4 outstat=factor;
            var CG:;
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
            merge factor_all_genes_communal (in=in1) flip (in=in2);
            by _name_;
            rename col1 = num_factor_&num;
            run;

        data factor_all_genes_communal;
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
    %iter_nfact(31);
    %iter_nfact(32);
    %iter_nfact(33);
    %iter_nfact(34);
    %iter_nfact(35);
    %iter_nfact(36);
    %iter_nfact(37);
    %iter_nfact(38);
    %iter_nfact(39);
    %iter_nfact(40);

    data SEM.dsx_factor_all_communal;
        set factor_all_genes_communal;
        run;


/* Run Factor Analysis Combined Genes */
    proc transpose data=SEM.dsrp_dsx_sbs_combine_sym(obs=1) out=flip;
        var CG:;
        run;
        
    data factor_combine_genes_communal;
        set flip;
        drop col1;
        run;
        
    proc sort data=factor_combine_genes_communal;
        by _name_;
        run;
    
    %macro iter_nfact(num);
        proc printto LOG="!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_combined_&num..log" PRINT="!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_combined_&num..lst" NEW;
            run;
        proc factor data=SEM.dsrp_dsx_sbs_combine_sym simple method=PRINCIPAL scree rotate=varimax nfact=&num round flag=.4 outstat=factor;
            var CG:;
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
            merge factor_combine_genes_communal (in=in1) flip (in=in2);
            by _name_;
            rename col1 = num_factor_&num;
            run;

        data factor_combine_genes_communal;
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
    %iter_nfact(31);
    %iter_nfact(32);
    %iter_nfact(33);
    %iter_nfact(34);
    %iter_nfact(35);
    %iter_nfact(36);
    %iter_nfact(37);
    %iter_nfact(38);
    %iter_nfact(39);
    %iter_nfact(40);

    data SEM.dsx_factor_combine_communal;
        set factor_combine_genes_communal;
        run;

/* Clean up */
proc datasets nolist;
    delete com;
    delete factor;
    delete factor_all_genes_communal;
    delete factor_combine_genes_communal;
    delete flip;
    delete merged;
    run; quit;
