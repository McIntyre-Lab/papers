/*******************************************************************************
* Filename: cegV_simulation_make_all_list.sas
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

/* Calculate mean and variance for all genes */
    proc sort data=SEM.cegsV_by_gene_stack  ;
        by sym;
        run;

    proc means data=SEM.cegsV_by_gene_stack  noprint;
        by sym;
        output out=all_summary mean(exp)=mean var(exp)=variance;
        run;

    ods graphics on;
    ods pdf file='!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_all_genes_distributions.pdf';
    title 'Distirubtion of Means; All';
    proc sgplot data=all_summary;
        density mean / type=kernel;
        run;
        
    title 'Distirubtion of Variances; All';
    proc sgplot data=all_summary;
        density variance / type=kernel;
        run;

    proc means data=all_summary min mean median max;
    var mean variance;
    run;

    ods pdf close;
    ods graphics off;

    data all_summary2;
        set all_summary;
        drop _type_ _freq_;
        run;

/* Export mean and variances */
    proc export data=all_summary2 outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_all_mean_and_variances.csv' dbms=csv replace;
        putnames=yes;
        run;
