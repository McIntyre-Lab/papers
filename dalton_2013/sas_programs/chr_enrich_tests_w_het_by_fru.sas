/********************************************************************************
* This script performs a fisher's exact test for x chromosome enrichment
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';

%macro enrichment_test(chrom,type,sex);

    proc freq data=FRU.flag_ind_rep_w_het_&sex.;
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
            length chrom $6.;
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
            length chrom $6.;
            set &sex._&chrom._&type._table;
            sex   = "&sex";
            name  = "&type";
            chrom = "&chrom";
            run;
    %end;

    proc datasets nolist;
        delete &sex._&chrom._&type._table;
        delete c&chrom._chisq_&type.;
        run;

%mend enrichment_test;

%enrichment_test(x     , a_ind   , male);
%enrichment_test(xhet  , a_ind   , male);
%enrichment_test(2L    , a_ind   , male);
%enrichment_test(2Lhet , a_ind   , male);
%enrichment_test(2R    , a_ind   , male);
%enrichment_test(2Rhet , a_ind   , male);
%enrichment_test(3L    , a_ind   , male);
%enrichment_test(3Lhet , a_ind   , male);
%enrichment_test(3R    , a_ind   , male);
%enrichment_test(3Rhet , a_ind   , male);
%enrichment_test(4     , a_ind   , male);
%enrichment_test(x     , a_rep , male);
%enrichment_test(xhet  , a_rep , male);
%enrichment_test(2L    , a_rep , male);
%enrichment_test(2Lhet , a_rep , male);
%enrichment_test(2R    , a_rep , male);
%enrichment_test(2Rhet , a_rep , male);
%enrichment_test(3L    , a_rep , male);
%enrichment_test(3Lhet , a_rep , male);
%enrichment_test(3R    , a_rep , male);
%enrichment_test(3Rhet , a_rep , male);
%enrichment_test(4     , a_rep , male);

%enrichment_test(x     , b_ind   , male);
%enrichment_test(xhet  , b_ind   , male);
%enrichment_test(2L    , b_ind   , male);
%enrichment_test(2Lhet , b_ind   , male);
%enrichment_test(2R    , b_ind   , male);
%enrichment_test(2Rhet , b_ind   , male);
%enrichment_test(3L    , b_ind   , male);
%enrichment_test(3Lhet , b_ind   , male);
%enrichment_test(3R    , b_ind   , male);
%enrichment_test(3Rhet , b_ind   , male);
%enrichment_test(4     , b_ind   , male);
%enrichment_test(x     , b_rep , male);
%enrichment_test(xhet  , b_rep , male);
%enrichment_test(2L    , b_rep , male);
%enrichment_test(2Lhet , b_rep , male);
%enrichment_test(2R    , b_rep , male);
%enrichment_test(2Rhet , b_rep , male);
%enrichment_test(3L    , b_rep , male);
%enrichment_test(3Lhet , b_rep , male);
%enrichment_test(3R    , b_rep , male);
%enrichment_test(3Rhet , b_rep , male);
%enrichment_test(4     , b_rep , male);

%enrichment_test(x     , c_ind   , male);
%enrichment_test(xhet  , c_ind   , male);
%enrichment_test(2L    , c_ind   , male);
%enrichment_test(2Lhet , c_ind   , male);
%enrichment_test(2R    , c_ind   , male);
%enrichment_test(2Rhet , c_ind   , male);
%enrichment_test(3L    , c_ind   , male);
%enrichment_test(3Lhet , c_ind   , male);
%enrichment_test(3R    , c_ind   , male);
%enrichment_test(3Rhet , c_ind   , male);
%enrichment_test(4     , c_ind   , male);
%enrichment_test(x     , c_rep , male);
%enrichment_test(xhet  , c_rep , male);
%enrichment_test(2L    , c_rep , male);
%enrichment_test(2Lhet , c_rep , male);
%enrichment_test(2R    , c_rep , male);
%enrichment_test(2Rhet , c_rep , male);
%enrichment_test(3L    , c_rep , male);
%enrichment_test(3Lhet , c_rep , male);
%enrichment_test(3R    , c_rep , male);
%enrichment_test(3Rhet , c_rep , male);
%enrichment_test(4     , c_rep , male);

%enrichment_test(x     , a_ind   , female);
%enrichment_test(xhet  , a_ind   , female);
%enrichment_test(2L    , a_ind   , female);
%enrichment_test(2Lhet , a_ind   , female);
%enrichment_test(2R    , a_ind   , female);
%enrichment_test(2Rhet , a_ind   , female);
%enrichment_test(3L    , a_ind   , female);
%enrichment_test(3Lhet , a_ind   , female);
%enrichment_test(3R    , a_ind   , female);
%enrichment_test(3Rhet , a_ind   , female);
%enrichment_test(4     , a_ind   , female);
%enrichment_test(x     , a_rep , female);
%enrichment_test(xhet  , a_rep , female);
%enrichment_test(2L    , a_rep , female);
%enrichment_test(2Lhet , a_rep , female);
%enrichment_test(2R    , a_rep , female);
%enrichment_test(2Rhet , a_rep , female);
%enrichment_test(3L    , a_rep , female);
%enrichment_test(3Lhet , a_rep , female);
%enrichment_test(3R    , a_rep , female);
%enrichment_test(3Rhet , a_rep , female);
%enrichment_test(4     , a_rep , female);

%enrichment_test(x     , b_ind   , female);
%enrichment_test(xhet  , b_ind   , female);
%enrichment_test(2L    , b_ind   , female);
%enrichment_test(2Lhet , b_ind   , female);
%enrichment_test(2R    , b_ind   , female);
%enrichment_test(2Rhet , b_ind   , female);
%enrichment_test(3L    , b_ind   , female);
%enrichment_test(3Lhet , b_ind   , female);
%enrichment_test(3R    , b_ind   , female);
%enrichment_test(3Rhet , b_ind   , female);
%enrichment_test(4     , b_ind   , female);
%enrichment_test(x     , b_rep , female);
%enrichment_test(xhet  , b_rep , female);
%enrichment_test(2L    , b_rep , female);
%enrichment_test(2Lhet , b_rep , female);
%enrichment_test(2R    , b_rep , female);
%enrichment_test(2Rhet , b_rep , female);
%enrichment_test(3L    , b_rep , female);
%enrichment_test(3Lhet , b_rep , female);
%enrichment_test(3R    , b_rep , female);
%enrichment_test(3Rhet , b_rep , female);
%enrichment_test(4     , b_rep , female);

%enrichment_test(x     , c_ind   , female);
%enrichment_test(xhet  , c_ind   , female);
%enrichment_test(2L    , c_ind   , female);
%enrichment_test(2Lhet , c_ind   , female);
%enrichment_test(2R    , c_ind   , female);
%enrichment_test(2Rhet , c_ind   , female);
%enrichment_test(3L    , c_ind   , female);
%enrichment_test(3Lhet , c_ind   , female);
%enrichment_test(3R    , c_ind   , female);
%enrichment_test(3Rhet , c_ind   , female);
%enrichment_test(4     , c_ind   , female);
%enrichment_test(x     , c_rep , female);
%enrichment_test(xhet  , c_rep , female);
%enrichment_test(2L    , c_rep , female);
%enrichment_test(2Lhet , c_rep , female);
%enrichment_test(2R    , c_rep , female);
%enrichment_test(2Rhet , c_rep , female);
%enrichment_test(3L    , c_rep , female);
%enrichment_test(3Lhet , c_rep , female);
%enrichment_test(3R    , c_rep , female);
%enrichment_test(3Rhet , c_rep , female);
%enrichment_test(4     , c_rep , female);

/* Export Chi-SQ Tables and Tests */

%macro export_chi(type);
    proc export data=chisq_&type
                outfile = "!MCLAB/arbeitman_fru_network/reports_external/chrom_enrichment_&type._w_het_by_fru.csv"
                dbms=CSV replace;
                putnames = yes;
                run;
%mend export_chi;

%export_chi(tests);
%export_chi(tables);
