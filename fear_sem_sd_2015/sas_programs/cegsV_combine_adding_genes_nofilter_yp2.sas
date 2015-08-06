/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegs_adding_genes_yp2/sas_data';

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

proc printto log='/home/jfear/tmp/cegv_yp2.log' new; run;
%iterdataset(dataset=SEM.cegsV_gene_list, function=%nrstr(%combine_bic(&symbol)));
proc printto log=log new; run;

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

data SEM.cegsV_ag_yp2_stack_bic;
    set combined;
    run;

proc transpose data=SEM.cegsV_ag_yp2_stack_bic out=flip;
    by gene;
    var BIC;
    id model;
    run;

data SEM.cegsV_ag_yp2_sbs_bic;
    retain gene Model_baseline Model_1 Model_2 Model_3 Model_4 Model_5 Model_6
    Model_7 Model_8 Model_9 Model_10 Model_11 Model_12 Model_13 Model_14
    Model_15 Model_16 Model_17 Model_18 Model_19 Model_20 Model_21 Model_22
    Model_23 Model_24 Model_25 Model_26 Model_27 Model_28 Model_29 Model_30
    Model_31 Model_32 Model_33 Model_34 Model_35 Model_36 Model_37
    ;
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
    delete gene_count;
    run; quit;
