/********************************************************************************
* This script performs a fisher's exact test for x chromosome enrichment
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';

%macro enrichment_test(chrom,type);

    proc freq data=FRU.flag_x_ind_rep_null_male;
        tables flag_&type.*flag_&chrom._chrom /chisq expected outexpect sparse out=c&chrom._&type._table;
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
            retain chrom name;
            set chisq_tests (in=in1)
                c&chrom._chisq_&type. (in=in2);
            if in2 then do;
                chrom = "&chrom";
                name = "&type";
            end;
            run;
    %end;
    %else %do;
        data chisq_tests;
            length name $9.;
            length chrom $5.;
            set c&chrom._chisq_&type.;
            chrom = "&chrom";
            name = "&type";
            run;
    %end;

    /* Format Chi-SQ Contingency Tables and Combine */
    data c&chrom._&type._table;
        set c&chrom._&type._table;
        rename flag_&type = flag_sig;
        rename flag_&chrom._chrom = flag_chrom;
        run;

    %if %sysfunc(exist(chisq_tables)) %then %do;
        data chisq_tables;
            retain name chrom;
            set chisq_tables (in=in1)
            c&chrom._&type._table (in=in2);
            if in2 then do;
                name  = "&type";
                chrom = "&chrom";
            end;
            run;
    %end;
    %else %do;
        data chisq_tables;
            length name $9.;
            length chrom $2.;
            set c&chrom._&type._table;
            name  = "&type";
            chrom = "&chrom";
            run;
    %end;

    proc datasets nolist;
        delete c&chrom._&type._table;
        delete c&chrom._chisq_&type.;
        run;quit;

%mend enrichment_test;

%enrichment_test(x,induced);
%enrichment_test(2L,induced);
%enrichment_test(2R,induced);
%enrichment_test(3L,induced);
%enrichment_test(3R,induced);
%enrichment_test(4,induced);
%enrichment_test(x,repressed);
%enrichment_test(2L,repressed);
%enrichment_test(2R,repressed);
%enrichment_test(3L,repressed);
%enrichment_test(3R,repressed);
%enrichment_test(4,repressed);

/* Export Chi-SQ Tables and Tests */

%macro export_chi(type);
    proc export data=chisq_&type
                outfile = "!MCLAB/arbeitman_fru_network/reports_external/chrom_enrichment_&type._null.csv"
                dbms=CSV replace;
                putnames = yes;
                run;
%mend export_chi;

%export_chi(tests);
%export_chi(tables);
