/*******************************************************************************
* Filename: cegV_data_prep_added_genes_for_mmc.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Output the 1390 genes for import into MMC. Genes are rows and
* lines are columns.
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

/* Transpose so lines are columns */
proc sort data=pos;
    by symbol;
    run;

proc transpose data=pos out=flip;
    by symbol;
    var exp;
    id line;
    run;

data clean;
    set flip;
    drop _name_;
    run;

/* Export for MMC */
proc export data=clean outfile='/home/jfear/mclab/cegs_sem_sd_paper/exported_data/cegsV_added_genes_for_mmc.csv' dbms=csv replace;
    putnames=yes;
    run;

