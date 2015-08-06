/*******************************************************************************
* Filename: cegV_simulation_make_null_list.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Create a list of genes that are not associated with the sex
* hierarchy.  Merge on the added genes and keep only those genes not added to
* the SD pathway.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge the genes added to the SD pathway to full gene list */
    proc sort data=SEM.cegsV_ag_yp2_added_genes;
        by gene;
        run;

    proc sort data=SEM.cegsV_gene_list;
        by symbol;
        run;

    data pos null;
        merge SEM.cegsV_ag_yp2_added_genes (in=in1 rename=(gene = symbol)) SEM.cegsV_gene_list (in=in2);
        by symbol;
        if in1 then output pos;
        else if in2 then output null;
        run;

/* Drop genes in the SD pathway from the null list */
    data null2;
        set null;
        if symbol = 'Spf45' then delete;
        if symbol = 'Sxl' then delete;
        if symbol = 'Yp2' then delete;
        if symbol = 'Yp3' then delete;
        if symbol = 'dsx' then delete;
        if symbol = 'fl_2_d' then delete;
        if symbol = 'fru' then delete;
        if symbol = 'her' then delete;
        if symbol = 'ix' then delete;
        if symbol = 'snf' then delete;
        if symbol = 'tra' then delete;
        if symbol = 'tra2' then delete;
        if symbol = 'vir' then delete;
        drop model;
        run;

/* Pull data for the NULL genes */
    proc sort data=SEM.cegsV_by_gene_stack;
        by sym;
        run;

    data merged2;
        merge null2 (in=in1) SEM.cegsV_by_gene_stack (in=in2 rename=(sym=symbol));
        by symbol;
        if in1 and in2;
        run;

/* Calculate mean and variance for each gene */
    proc means data=merged2 noprint;
        by symbol;
        output out=summary mean(exp)=mean var(exp)=variance;
        run;

    ods graphics on;
    ods pdf file='!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_null_distributions.pdf';
    title 'Distirubtion of Means; Null';
    proc sgplot data=summary;
        density mean / type=kernel;
        run;
        
    title 'Distirubtion of Variances; Null';
    proc sgplot data=summary;
        density variance / type=kernel;
        run;

    proc means data=summary min mean median max;
    var mean variance;
    run;

    ods pdf close;
    ods graphics off;

    data summary2;
        set summary;
        drop _type_ _freq_;
        run;

/* Export mean and variances */
    proc export data=summary2 outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_null_mean_and_variances.csv' dbms=csv replace;
        putnames=yes;
        run;
