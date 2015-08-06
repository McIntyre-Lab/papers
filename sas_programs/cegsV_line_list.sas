/*******************************************************************************
* Filename: cegsV_line_list.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Create a list of lines that we are using in the cegs SD paper
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data SEM.cegsV_line_list;
    set SEM.cegsV_by_gene_sbs;
    keep line;
    run;

proc export data=SEM.cegsV_line_list outfile='!MCLAB/cegs_sem_sd_paper/design_file/cegsV_line_list.csv' dbms=csv replace;
    putnames=no;
    run;
