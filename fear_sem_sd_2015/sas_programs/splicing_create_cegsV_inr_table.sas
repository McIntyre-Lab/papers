/*******************************************************************************
* Filename: misc_test_create_cegsV_inr_table.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Build a simple sbs table with InR fusions and their normalized
* values. Note that several InR fusions are missing because they were removed
* during the normalization process.
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Get the InR gene region */
    data _null_;
        set DMEL551.symbol2coord;
        where symbol eq 'InR';
        call symput('chrom', chrom);
        call symput('start', start);
        call symput('end', end);
        run;

/* Pull fusions within that region */
    data fusList;
        set DMEL551.fb551_si_fusions_unique_flagged;
        where chrom eq "&chrom" and start ge &start and end le &end;
        keep fusion_id;
        run; *18 obs including CR43653;

/* Pull out above fusions from data */
    proc sort data=fusList;
        by fusion_id;
        run;

    proc sort data=SEM.cegs_virgin_norm_cent;
        by fusion_id;
        run;

    data fusTable;
        merge fusList (in=in1) SEM.cegs_virgin_norm_cent (in=in2);
        by fusion_id;
        if in1;
        keep fusion_id line mean_exp;
        run;

/* Create table with cols=fusion_id, row=line, value=mean_exp */
    proc sort data=fusTable;
        by line fusion_id;
        run;

    proc transpose data=fusTable out=flip;
        by line;
        var mean_exp;
        id fusion_id;
        run;

    data SEM.inr_coverage;
        retain line F57608_SI F57609_SI S57592_SI S57593_SI S57594_SI S57595_SI
        S57596_SI S57597_SI S57598_SI S57599_SI S57600_SI S57601_SI S57602_SI
        S57603_SI S57604_SI S57605_SI S57606_SI S57607_SI
        ;
        set flip;
        if line eq ' ' then delete;
        drop _name_;
        run;

/* Clean up */
proc datasets ;
    delete flip;
    delete fuslist;
    delete fustable;
    run; quit;

