/*******************************************************************************
* Filename: enrichment_go_data_prep.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I will do go enrichment in JMP Genomics, but need to prep the
* dataset here.
*
*******************************************************************************/

/* Libraries
    libname SEM '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Merge data to go terms by fbgn */
proc sort data=SEM.cegsV_ag_w_flags_bic12;
    by primary_fbgn;
    run;

proc sort data=DMEL551.genes2go_nodups;
    by fbgn;
    run;

data merged;
    merge SEM.cegsV_ag_w_flags_bic12 (in=in1) DMEL551.genes2go_nodups (in=in2 rename=(fbgn=primary_fbgn));
    by primary_fbgn;
    if in1;
    run;

data SEM.cegsV_ag_w_flags_bic12_go;
    set merged;
    run;
