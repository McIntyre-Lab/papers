/*******************************************************************************
* Filename: enrichmentTest.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: A macro to do enrichment analysis.
*
* INPUTS:
*   bg = Dataset containg your category information (i.e., gene list merged to
*        go terms)
*
*
*   sigBinFlag = Dataset containg your binary flag of significance
*
*   bgID = name of primary key for merging to flags (i.e., fbgn)
*
*   flagID = name of primary key for merging to bg file (i.e., fbgn)
*
*   category = name of the column with categories in background dataset
*
*   delimiter = delimiter if there are multiple categories in the category column
*
* OUTPUTS: 
*   out = name of the output dataset containing 
*
*******************************************************************************/

%macro enrichmentTest(bg, sigBinFlag, bgID, flagID, category, delimiter, out);
    /*
    * bg = Dataset containg your category information (i.e., gene list merged to go terms)
    * sigBinFlag = Dataset containg your binary flag of significance
    *
    * bgID = name of primary key for merging to flags (i.e., fbgn)
    * flagID = name of primary key for merging to bg file (i.e., fbgn)
    *
    * category = name of the column with categories in background dataset
    * delimiter = delimiter if there are multiple categories in the category column
    *
    * out = name of the output dataset
    */

    /* split concatentated terms and make a 0|1 flag */
        proc sort data=&bg;
            by &bgID;
            run;

        data term;
            set &bg;
            by &bgID;
            do i=1 by 1 while(scan(&category, i, &delimiter)^='');
            term = scan(&category, i, &delimiter);
            flag_term = 1;
            output;
            end;
            rename &bgID = ID;
            keep &bgID term flag_term;
            run;

        %* create a per gene flag for each go term;
        proc transpose data=term out=flip;
            by ID;
            var flag_term;
            id term;
            run;

        proc transpose data=flip out=flip_back;
            by ID;
            var :;
            id _name_;
            run;

        data background;
            set flip_back;
            if flag_term = ' ' then flag_term = 0;
            term = tranwrd(_name_, '_', ':');
            drop _name_;
            run;

    /* Merge on flags */
        proc sort data=&sigBinFlag;
            by &flagID;
            run;

        proc sort data=background;
            by ID;
            run;

        data merged;
            merge &sigBinFlag (in=in1 rename=(&flagID = ID)) background (in=in2);
            by ID;
            if in1 then flag = 1;
            else flag = 0;
            if in1 and not in2 then flag_term = 0;
            run;

    /* Enrichment Analsyis */
        proc sort data=merged;
            by term;
            run;

        proc freq data=merged noprint;
            by term;
            tables flag*flag_term /fisher out=freqs;
            exact fisher;
            output out=tests fisher;
            run;

        data enrich;
            set tests;
            label XP2_FISH = ' ';
            keep term XP2_FISH;
            run;

    /* FDR correction */
        data raw_p;
            set enrich;
            rename XP2_FISH = raw_p;
            run;

        proc multtest inpvalues=raw_p fdr out=fdr noprint;
            run; quit;

    /* Summarize Counts */
        proc sort data=freqs;
            by term;
            run;

        data cat;
            set freqs;
            if term eq '' then delete;
            category = 'class_' || strip(flag) || '_' || strip(flag_term);
            run;

        proc transpose data=cat out=catFlip;
            by term;
            var count;
            id category;
            run;

        data catTotal;
            set catFlip;
            if category_0_0 = '' then category_0_0 = 0;
            if category_1_0 = '' then category_1_0 = 0;
            if category_0_1 = '' then category_0_1 = 0;
            if category_1_1 = '' then category_1_1 = 0;
            class_0_total = sum(category_0_0, category_0_1);
            class_0_cnt = category_0_1;
            class_0_prop = class_0_cnt / class_0_total;
            class_1_total = sum(category_1_0, category_1_1);
            class_1_cnt = category_1_1;
            class_1_prop = class_1_cnt / class_1_total;
            keep term class_:;
            run;

    /* Merge and output */
        data enrichment_results;
            merge fdr (in=in1) catTotal (in=in2);
            by term;
            label fdr_p = ' ';
            run;

        proc sort data=enrichment_results;
            by fdr_p;
            run;

    /* Clean up */
        proc datasets ;
            delete term;
            delete flip;
            delete flip_back;
            delete background;
            delete merged;
            delete freqs;
            delete tests;
            delete enrich;
            delete raw_p;
            delete fdr;
            delete cat;
            delete catflip;
            delete cattotal;
            run; quit;
%mend;
