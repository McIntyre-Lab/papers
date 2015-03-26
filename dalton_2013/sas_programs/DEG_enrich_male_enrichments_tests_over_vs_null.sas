%macro loop(comp1,comp2);
    proc freq data=FRU.flag_male_status;
        table flag_&comp1*flag_&comp2 /chisq expected outexpect sparse out=fru_&comp1._vs_&comp2._table;
        exact fisher chisq;
        output out=fru_&comp1._vs_&comp2._test fisher chisq;
        run;

    /* Format Chi-sq tests and combine */
        data fru_&comp1._vs_&comp2._test2;
            set fru_&comp1._vs_&comp2._test;
            label XP2_FISH = 'Fisher raw P-value (2-tail)';
            keep N _PCHI_ DF_PCHI XP_PCHI XP2_FISH;
            run;

        %if %sysfunc(exist(chisq_tests)) %then %do;
            data chisq_tests;
                retain comp1 comp2;
                set chisq_tests (in=in1)
                    fru_&comp1._vs_&comp2._test2 (in=in2);
                if in2 then do;
                    comp1 = "&comp1";
                    comp2 = "&comp2";
                end;
                run;
        %end;
        %else %do;
            data chisq_tests;
                length comp1 $4.;
                length comp2 $4.;
                set fru_&comp1._vs_&comp2._test2;
                comp1 = "&comp1";
                comp2 = "&comp2";
                run;
        %end;

    /* Format Chi-SQ Contingency Tables and Combine */
        data fru_&comp1._vs_&comp2._table2;
            set fru_&comp1._vs_&comp2._table;
            rename flag_&comp1 = flag_comp1;
            rename flag_&comp2 = flag_comp2;
            run;

        %if %sysfunc(exist(chisq_tables)) %then %do;
            data chisq_tables;
                retain comp1 comp2;
                set chisq_tables (in=in1)
                    fru_&comp1._vs_&comp2._table2 (in=in2);
                if in2 then do;
                    comp1  = "&comp1";
                    comp2 = "&comp2";
                end;
                run;
        %end;
        %else %do;
            data chisq_tables;
                length comp1 $4.;
                length comp2 $4.;
                set fru_&comp1._vs_&comp2._table2;
                comp1  = "&comp1";
                comp2 = "&comp2";
                run;
        %end;

    /* Clean up */
    proc datasets nolist;
        delete fru_&comp1._vs_&comp2._test fru_&comp1._vs_&comp2._test2;
        delete data fru_&comp1._vs_&comp2._table fru_&comp1._vs_&comp2._table2; 
        run;quit;
%mend;
%loop(a,b);
%loop(a,c);
%loop(a,null);
%loop(b,c);
%loop(b,null);
%loop(c,null);


/* Export Chi-SQ Tables and Tests */

%macro export_chi(type);
    proc export data=chisq_&type
                outfile = "!MCLAB/arbeitman/arbeitman_fru_network/reports_external/male_ind_rep_enrichment_&type..csv"
                dbms=CSV replace;
                putnames = yes;
                run;
%mend export_chi;

%export_chi(tests);
%export_chi(tables);
