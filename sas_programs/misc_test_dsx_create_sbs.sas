/********************************************************************************
* Using the DSX datasets, I have created a gene list with genes repressed
* when dsxF is removed (dsxNull females).
********************************************************************************/
/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
    libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Merge on Fb r5.48 gene annotations */
    * The dsx analysis was done using FB r5.30. In addition the dsx gene lists
    * contain only FBgn number. In order to merge to the dsrp dataset I need to
    * merge on gene annotations (FB r5.48) and use CGnumber.
    ;

    /** Merge on Fb5.30 symbol **/
        * Sometimes the primary FBgn number changes between flybase relases, so
        * merge on gene symbol
        ;
        data dsxNull;
            set SEM.dsxNullf_repressed;
            rename fbgn_cat = primary_fbgn;
            run;

        proc sort data=dsxNull;
            by primary_fbgn;
            run;

        proc sort data=dmel530.symbol2cg;
            by primary_fbgn;
            run;

        data merged oops;
            merge dmel530.symbol2cg (in=in1) dsxNull (in=in2);
            by primary_fbgn;
            if in1 and in2 then output merged;
            else if in2 then output oops;
            run;

        data merged_rename;
            set merged;
            rename annotation_id = cgnumber;
            * Some of the genes would not line up with the below merge. So I
            * went to the flybase archive and looked up the ones that did not
            * merge properly and fixed them here.;
            if symbol eq 'bru-3' then annotation_id = 'CG43744';
            if symbol eq 'cac' then annotation_id = 'CG43368';
            if symbol eq 'CG15928' then delete; * this gene was merged by flybase with cac;
            if symbol eq 'btsz' then annotation_id = 'CG44012';
            if symbol eq 'CG34400' then annotation_id = 'CG43749';
            if symbol eq 'CG3563' then annotation_id = 'CG43162';
            if symbol eq 'CG42737' then delete; * Fb removed this gene;
            if symbol eq 'cngl' then annotation_id = 'CG43395';
            run; * Now I only have 206 obs instead of 208 because I delete two.;

    /** Merge Fb5.48 isoform information to symbol **/
        * The dsrp dataset is based on the protein isoform cgnumbers. I need to
        * merge on isoform information prior to merging to dsrp.
        ;
        proc sort data= merged_rename;
            by cgnumber;
            run;

        proc sort data=dmel548.symbol2proteincg;
            by cgnumber;
            run;

        data merged2 oops2;
            merge dmel548.symbol2proteincg (in=in1) merged_rename (in=in2);
            by cgnumber;
            if in1 and in2 then output merged2;
            else if in2 then output oops2;
            run;

        /* Test that I get my 206 back 

            proc sort data=merged2 nodupkey out=nodup;
            by cgnumber;
            run;

            Success 206!!!
        */

        data merged2_rename;
            set merged2;
            rename protein_cg_number = cgnumber;
            keep protein_cg_number;
            run;

        /* Another check 

            data test;
                set merged2_rename;
                num = index(cgnumber, '_');
                sym = substr(cgnumber, 1, num-1);
                keep sym;
                run;

            proc sort data=test nodupkey;
                by sym;
                run;

            Success 206!!!

        */


    /** Merge Fb r5.48 gene list to dsrp dataset **/
        proc sort data=merged2_rename;
            by cgnumber;
            run;

        proc sort data=sem.dsrp_stack;
            by cgnumber;
            run;

        data merged3 nodsrp;
            merge sem.dsrp_stack (in=in1) merged2_rename (in=in2);
            by cgnumber;
            if in1 and in2 then output merged3;
            else if in2 then output nodsrp;
            run;

    /* check how many genes are left in the DSRP dataset
            data test;
                set merged3;
                num = index(cgnumber, '_');
                sym = substr(cgnumber, 1, num-1);
                keep sym;
                run;

            proc sort data=test nodupkey;
                by sym;
                run;

            Only 35 Genes left :(

    */

/* transpose dataset into side-by-side and make perminant */
    proc sort data=merged3;
        by patRIL matRIL;
        run;

    proc transpose data=merged3 out=flip;
        by patRIL matRIL;
        id cgnumber;
        var col1;
        run;

    data SEM.dsrp_dsx_sbs_symbol;
        set flip;
        drop _name_;
        run;

/* Create cgisoform list */
    data SEM.dsxNullf_repressed_cgnumber;
        set merged3;
        keep cgnumber;
        run;

/* Clean up */
    proc datasets nolist;
        delete dsxnull;
        delete flip;
        delete merged;
        delete merged2;
        delete merged2_rename;
        delete merged3;
        delete merged3_rename;
        delete merged_rename;
        delete oops;
        delete oops2;
        delete nodsrp;
        delete nodup;
        delete test;
    run; quit;
