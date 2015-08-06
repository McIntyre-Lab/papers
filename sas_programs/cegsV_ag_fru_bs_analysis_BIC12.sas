/*******************************************************************************
* Filename: cegsV_ag_fru_bs_analysis_BIC12.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Using the previous Dalton 2013 Fru BS data, are the 9 genes that
* are best downstream of fru have the fru BS?
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Grab genes added best downstream of fru */
    data dsfru;
        set SEM.cegsV_ag_yp2_flag_ds_fru_BIC12;
        where flag_best_fru = 1;
        *where flag_all_fru = 1;
        drop primary_fbgn;
        run;

/* Convert to Fb530 */
    * Motifs were called on Fb530, so I need to convert the primary fbgns;

    proc sort data=dsfru;
        by symbol;
        run;

    proc sort data=DMEL530.symbol2fbgn;
        by symbol;
        run;

    data dsfru530;
        merge dsfru (in=in1) DMEL530.symbol2fbgn (in=in2);
        by symbol;
        if in1;
        run;

    proc sort data=dsfru530;
        by primary_fbgn;
        run;


/* Grab Motifs from Fru network project */
    data motif;
        set FRU.motif_flags_and_cnts;
        keep primary_fbgn flag_fru_a_motif flag_fru_b_motif flag_fru_c_motif;
        run;

    proc sort data=motif;
        by primary_fbgn;
        run;

    data merged;
        merge motif (in=in1) dsfru530 (in=in2);
        by primary_fbgn;
        if in2;
        flag_fru_bs = flag_fru_a_motif + flag_fru_b_motif + flag_fru_c_motif;
        drop flag_fru_a_motif flag_fru_b_motif flag_fru_c_motif flag_all_fru;
        run;

/* 
data tmp;
    set merged;
    if flag_fru_bs gt 0;
    run;
*/

/* Export list with motif count */
    proc export data=merged outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_ds_fru_binding_site_BIC12.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
proc datasets ;
    delete dsfru;
    delete dsfru530;
    delete merged;
    delete motif;
    run; quit;

