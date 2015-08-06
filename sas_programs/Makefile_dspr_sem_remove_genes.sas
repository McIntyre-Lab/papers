/*******************************************************************************
* Filename: Makefile_dspr_sem_remove_genes.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Remove each gene in the sex hierarchy one at a time and see how
* that affects BIC.
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Combine Remove Gene Simulation */
    * I ran 100 simulations on the HPC: (1) take original data and no
    * covariance model and simulate a new dataset, (2) remove each gene from
    * the model one at a time. Now I need to combine these results into a single table.
    *
    * INPUT: !HOME/tmp/dspr_removing_genes_simulation/&i/fitstat
    * 
    * DATASET:
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_sem_remove_genes_combine_simulation.sas';

