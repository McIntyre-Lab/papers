/********************************************************************************
* This script performs a fisher's exact test for x chromosome enrichment
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';

%macro enrichment_test(chrom,type,sex);

    proc freq data=FRU.flag_x_induced_repressed_&sex.;
        tables flag_&type.*flag_&chrom._chrom /chisq expected outexpect sparse out=&sex._&chrom._&type._table;
        exact fisher chisq;
        output out=c&chrom._chisq_&type. fisher chisq;
        run;
        
    /* Format Chi-SQ tests and combine */
    data c&chrom._chisq_&type.;
        set c&chrom._chisq_&type.;
        label XP2_FISH = 'Fisher raw P-value (2-tail)';
        keep N _PCHI_ DF_PCHI XP_PCHI XP2_FISH ;
        run;


    %if %sysfunc(exist(chisq_tests)) %then %do;
        data chisq_tests;
            retain chrom name sex;
            set chisq_tests (in=in1)
                c&chrom._chisq_&type. (in=in2);
            if in2 then do;
                chrom = "&chrom";
                name = "&type";
                sex = "&sex";
            end;
            run;
    %end;
    %else %do;
        data chisq_tests;
            length name $9.;
            length sex $6.;
            set c&chrom._chisq_&type.;
            chrom = "&chrom";
            name = "&type";
            sex = "&sex";
            run;
    %end;

    /* Format Chi-SQ Contingency Tables and Combine */
    data &sex._&chrom._&type._table;
        set &sex._&chrom._&type._table;
        rename flag_&type = flag_sig;
        rename flag_&chrom._chrom = flag_chrom;
        run;

    %if %sysfunc(exist(chisq_tables)) %then %do;
        data chisq_tables;
            retain sex name chrom;
            set chisq_tables (in=in1)
            &sex._&chrom._&type._table (in=in2);
            if in2 then do;
                sex   = "&sex";
                name  = "&type";
                chrom = "&chrom";
            end;
            run;
    %end;
    %else %do;
        data chisq_tables;
            length name $9.;
            length sex $6.;
            length chrom $2.;
            set &sex._&chrom._&type._table;
            sex   = "&sex";
            name  = "&type";
            chrom = "&chrom";
            run;
    %end;

    proc datasets;
        delete &sex._&chrom._&type._table;
        delete c&chrom._chisq_&type.;
        run;

%mend enrichment_test;

%enrichment_test(x,induced,male);
%enrichment_test(2L,induced,male);
%enrichment_test(2R,induced,male);
%enrichment_test(3L,induced,male);
%enrichment_test(3R,induced,male);
%enrichment_test(4,induced,male);
%enrichment_test(x,repressed,male);
%enrichment_test(2L,repressed,male);
%enrichment_test(2R,repressed,male);
%enrichment_test(3L,repressed,male);
%enrichment_test(3R,repressed,male);
%enrichment_test(4,repressed,male);

%enrichment_test(x,induced,female);
%enrichment_test(2L,induced,female);
%enrichment_test(2R,induced,female);
%enrichment_test(3L,induced,female);
%enrichment_test(3R,induced,female);
%enrichment_test(4,induced,female);
%enrichment_test(x,repressed,female);
%enrichment_test(2L,repressed,female);
%enrichment_test(2R,repressed,female);
%enrichment_test(3L,repressed,female);
%enrichment_test(3R,repressed,female);
%enrichment_test(4,repressed,female);

/* Export Chi-SQ Tables and Tests */

%macro export_chi(type);
    proc export data=chisq_&type
                outfile = "!MCLAB/arbeitman_fru_network/reports_external/chrom_enrichment_&type..csv"
                dbms=CSV replace;
                putnames = yes;
                run;
%mend export_chi;

%export_chi(tests);
%export_chi(tables);
