/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

*libname addgen '!MCLAB/cegs_sem_sd_paper/adding_genes/newlink/sas_data';
libname addgen '/home/jfear/tmp/newlink/sas_data';

* Sex det genes;
    data sexdet;
        input symbol $;
        cards; 
            sxl
            tra
            fru
            Yp2
            fl_2_d
            Spf45
            snf
            vir
            tra2
            her
            ix
        ;
        run;

*Combine all genes into a single dataset;
%macro combine_bic(gene);
    %if %sysfunc(exist(addgen.&gene)) %then %do;
            
        proc append data=addgen.&gene base=combined;
            run;

    %end;
    %else %do;
        /* I did not run the key sex det genes because they were already in the
         * model. So they should not have a dataset. Check to make sure I got
         * everything elese. */
        data misgene;
            format gene $12.;
            gene = "&gene";
            run;

        proc append data=misgene base=checks;
            run;
    %end;
%mend;

%iterdataset(dataset=sexdet, function=%nrstr(%combine_bic(&symbol)));

* checks looks good;

/* The datasets in addgen do not contain a single gene. In other words I had
 * some appending issues. So now genes are in the datset multiple times. Just
 * sort nodupkey to get things right.
 */

/* Check for dups
proc sort data=combined nodupkey out=combined_uniq;
    by gene model;
    run;
*/

/* Create Perminant stacked datset */
proc sort data=combined;
    by gene;
    run;

proc freq data=combined noprint;
    table gene / out=gene_count;
    run;

data SEM.dsrp_add_newlink_stack_bic;
    set combined;
    run;

proc transpose data=combined out=flip;
    by gene;
    var BIC;
    id model;
    run;

data SEM.dsrp_add_newlink_sbs_bic;
    retain gene baseline Model_1 Model_2 Model_3;
    set flip;
    drop _name_;
    run;

/* Clean up */
proc datasets nolist;
    delete checks;
    delete combined;
    delete combined_uniq;
    delete flip;
    delete misgene;
    delete sbs;
    delete baseline;
    delete gene_count;
    delete sexdet;
    run; quit;
