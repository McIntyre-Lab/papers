/*******************************************************************************
* Filename: enrichment_ag_chromosome_enrichment.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Are genes added to the network enriched on the different
* chromosomes?
*
*******************************************************************************/

/* Libraries
    libname SEM '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/
/* Merge on Chromosome Information */
proc sort data=SEM.cegsV_ag_w_flags_bic12;
    by primary_fbgn;
    run;

proc sort data=DMEL551.fbgn2coord;
    by primary_fbgn;
    run;

data merged oops;
    merge SEM.cegsV_ag_w_flags_bic12 (in=in1) DMEL551.fbgn2coord (in=in2);
    by primary_fbgn;
    if in1 and in2 then output merged;
    else if in1 then output oops;
    drop start end strand ;
    run;

proc sort data=merged;
    by chrom;
    run;

data flag_chr;
    set merged;
    if chrom eq '2R' or chrom eq '2L' or chrom eq '2RHet' or chrom eq '2LHet' then flag_2 = 1; else flag_2 = 0;
    if chrom eq '3R' or chrom eq '3L' or chrom eq '3RHet' or chrom eq '3LHet' then flag_3 = 1; else flag_3 = 0;
    if chrom eq '4' then flag_4 = 1; else flag_4 = 0;
    if chrom eq 'X' or chrom eq 'XHet' then flag_x = 1; else flag_x = 0;
    run;


/* Chromosomal enrichment */
    %macro enrichment_test(chrom);

        proc freq data=flag_chr noprint;
            tables flag_expanded*flag_&chrom. /chisq expected outexpect sparse out=_&chrom._table;
            exact fisher chisq;
            output out=c&chrom._chisq fisher chisq;
            run;
            
        /* Format Chi-SQ tests and combine */
        data c&chrom._chisq;
            set c&chrom._chisq;
            label XP2_FISH = 'Fisher raw P-value (2-tail)';
            keep N _PCHI_ DF_PCHI XP_PCHI XP2_FISH ;
            run;


        %if %sysfunc(exist(chisq_tests)) %then %do;
            data chisq_tests;
                retain chrom ;
                set chisq_tests (in=in1)
                    c&chrom._chisq (in=in2);
                if in2 then do;
                    chrom = "&chrom";
                end;
                run;
        %end;
        %else %do;
            data chisq_tests;
                length chrom $2.;
                set c&chrom._chisq;
                chrom = "&chrom";
                run;
        %end;

        /* Format Chi-SQ Contingency Tables and Combine */
        data _&chrom._table;
            set _&chrom._table;
            rename flag_expanded = flag_sig;
            rename flag_&chrom. = flag_chrom;
            run;

        %if %sysfunc(exist(chisq_tables)) %then %do;
            data chisq_tables;
                retain chrom;
                set chisq_tables (in=in1)
                _&chrom._table (in=in2);
                if in2 then do;
                    chrom = "&chrom";
                end;
                run;
        %end;
        %else %do;
            data chisq_tables;
                length chrom $2.;
                set _&chrom._table;
                chrom = "&chrom";
                run;
        %end;

        proc datasets;
            delete _&chrom._table;
            delete c&chrom._chisq;
            run;

    %mend enrichment_test;

    %enrichment_test(2);
    %enrichment_test(3);
    %enrichment_test(4);
    %enrichment_test(x);

/* Export */
    proc export data=chisq_tests outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/enrichment/chromosome_enrichment_tests.csv' label dbms=csv replace ;
        putnames=yes;
        run;

    proc export data=chisq_tables outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/enrichment/chromosome_enrichment_tables.csv' label dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
proc datasets nolist;
    delete chisq_tables;
    delete chisq_tests;
    delete flag_chr;
    delete merged;
    delete oops;
    run; quit;

