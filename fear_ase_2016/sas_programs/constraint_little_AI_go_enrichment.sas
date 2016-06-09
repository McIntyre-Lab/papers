/*******************************************************************************
* Filename: constraint_little_AI_go_enrichment.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Perform a GO enrichment analysis on genes with little evidence
* for AI.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Create gene list of all genes that have at least 25 lines */
    data mated;
        set CEGS.mated_flag_ase_line_cnts;
        where total > 25;
        keep symbol_cat;
        run; * 4610;

    data virgin;
        set CEGS.virgin_flag_ase_line_cnts;
        where total > 25;
        keep symbol_cat;
        run; * 4844;

    data genes;
        set mated virgin;
        by symbol_cat;
        symbol = scan(symbol_cat, 1, '|');
        output;
        symbol = scan(symbol_cat, 2, '|');
        if symbol ne '' then output;
        drop symbol_cat;
        run;

    proc sort data=genes nodupkey;
        by symbol;
        run; *5150;

/* Merge on GO */
    * merge on FBGN;
    proc sort data=DMEL551.symbol2fbgn;
        by symbol;
        run;

    data wFBGN;
        merge genes (in=in1) DMEL551.symbol2fbgn (in=in2);
        by symbol;
        if in1 and in2;
        run;

    * merge on GO;
    proc sort data=wFBGN nodupkey;
        by primary_fbgn;
        run; *ok, 5150 obs no dups;

    proc sort data=DMEL551.genes2go_nodups;
        by fbgn;
        run;

    data go;
        retain fbgn symbol GO_number_biopro_cat GO_number_celcomp_cat
        GO_number_molfunc_cat GO_biological_process_cat
        GO_cellular_component_cat GO_molecular_function_cat;
        merge wFBGN (in=in1 rename=(primary_fbgn = fbgn)) DMEL551.genes2go_nodups (in=in2);
        by fbgn;
        if in1;
        run;

/* Create per gene flags */
    * Split concatenated symbols;
    data flags;
        set CEGS.flag_no_ai_mated_and_virgin;
        by symbol_cat;
        symbol = scan(symbol_cat, 1, '|');
        output;
        symbol = scan(symbol_cat, 2, '|');
        if symbol ne '' then output;
        drop symbol_cat;
        run;

    * merge on fbgn;
    proc sort data=flags nodupkey;
        by symbol;
        run; *319;

    proc sort data=DMEL551.symbol2fbgn;
        by symbol;
        run;

    data flagsFbgn;
        merge flags (in=in1) DMEL551.symbol2fbgn (in=in2);
        by symbol;
        if in1;
        run; *319;

%enrichmentTest(go, flagsFbgn, fbgn, primary_fbgn, go_number_biopro_cat, '|', biopro);

/* Clean up */
    proc datasets ;
        delete flags;
        delete flagsfbgn;
        delete genes;
        delete go;
        delete mated;
        delete out;
        delete virgin;
        delete wfbgn;
        run; quit;
