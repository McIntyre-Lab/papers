/*******************************************************************************
* Filename: cegsV_sem_remove_genes_combine_simulation.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I ran 100 simulations on the HPC: (1) take original data and no
* covariance model and simulate a new dataset, (2) remove each gene from the
* model one at a time. Now I need to combine these results into a single table.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Copy all of the data into a single data set */
%macro import_all_data();
    %do i=1 %to 100;
        libname curr "!HOME/tmp/cegs_removing_genes_simulation/&i";
        data curr;
            format sim best12.;
            set CURR.fitstat;
            sim = &i;
            run;

        proc append data=curr base=combined;
        run;
    %end;
%mend;
%import_all_data();

data better;
    set combined;
    by sim;
    flag + 0 ;
    if first.sim then flag=0;
    if gene_dropped eq 'baseline' then flag = 1;
    if flag eq 1 and gene_dropped ne 'baseline' then output;
    drop flag;
    run;

proc freq data=better;
table gene_dropped;
run;

/* Clean up */
proc datasets ;
    delete BETTER;
    delete COMBINED;
    delete CURR;
    run; quit;

