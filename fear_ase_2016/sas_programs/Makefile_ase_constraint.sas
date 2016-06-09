/*******************************************************************************
* Filename: Makefile_ase_constraint.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Looking at the sex determination pathway, we see that Sxl
* appears to be constrained to not have any AI. We are interested if this is an
* evolutionary process. Here I explore this hypothesis and try to decide if it
* is worth pursuing further.
*
*******************************************************************************/

/* Libraries */
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Enrichment Tests */
    * I have generated a set of gene lists that either have little AI or a lot
    * of AI. I want to see if these gene lists are enriched for anything
    * interesting.
    ;

    /* Little AI GO Enrichment */
        * Perform a GO enrichment analysis on genes with little evidence for
        * AI.
        *
        * INPUT: CEGS.mated_flag_ase_line_cnts
        *        CEGS.virgin_flag_ase_line_cnts
        *        DMEL551.symbol2fbgn
        *        DMEL551.genes2go_nodups
        *
        * DATASET: 
        *
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/constraint_little_AI_go_enrichment.sas';
