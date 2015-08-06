/*******************************************************************************
* Filename: cegV_simulation_added_distribution.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Create a list of genes that are not associated with the sex
* hierarchy.  Merge on the added genes and keep only those genes not added to
* the SD pathway.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge the genes added to the SD pathway to full gene list */
proc sort data=SEM.cegsV_ag_yp2_added_genes;
    by gene;
    run;

proc sort data=SEM.cegsV_by_gene_stack;
    by sym;
    run;

data pos;
    merge SEM.cegsV_ag_yp2_added_genes (in=in1 rename=(gene = symbol)) SEM.cegsV_by_gene_stack (in=in2 rename=(sym = symbol));
    by symbol;
    if in1 ;
    drop model;
    run;

/* Calculate mean and variance for each gene */
proc means data=pos noprint;
    by symbol;
    output out=summary mean(exp)=mean var(exp)=variance;
    run;

ods graphics on;
ods pdf file='!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_added_distributions.pdf';
title 'Distirubtion of Means; Added';
proc sgplot data=summary;
    density mean / type=kernel;
    run;
    
title 'Distirubtion of Variances; Added';
proc sgplot data=summary;
    density variance / type=kernel;
    run;

proc means data=summary min mean median max;
var mean variance;
run;

ods pdf close;
ods graphics off;

