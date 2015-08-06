/*******************************************************************************
* Filename: cegsV_random_subset.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: I want to run the CEGS data using half of the genotypes. If our
* hypothesis about the DSPR is true then we should not add any genes when using
* half the number of genotypes.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
*/

/* Take a random subset of cegs lines */
proc surveyselect data=SEM.cegsV_line_list method=SRS rep=1 sampsize=35
seed=12345 out=small;
id line;
run;

data SEM.cegsV_random_subset;
    set small;
    keep line;
    run;

/* Create sbs for running models */
proc sort data=SEM.cegsV_by_gene_sbs;
    by line;
    run;

proc sort data=SEM.cegsV_random_subset;
    by line;
    run;

data SEM.cegsV_by_gene_sbs_subset;
    merge SEM.cegsV_by_gene_sbs (in=in1) SEM.cegsV_random_subset (in=in2);
    by line;
    if in2;
    run;
