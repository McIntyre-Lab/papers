/********************************************************************************
* Create a gene list of genes with multiple sex determination genes in their
* primary neighborhood.
*
* libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
* libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
* filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
* options SASAUTOS=(sasautos mymacros);
********************************************************************************/

/* Import Neighborhood table */
    proc import datafile='!MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_neighborhood_analysis_v3.csv' dbms=csv out=neigh replace;
        getnames=yes;
        guessingrows=11064;
        run;

/* Keep only genes that have two or more sex determination genes in their primary neighborhood */
    data neigh2;
        retain FG cgnumber;
        set neigh;
        where primary_sex_det_genes ? '|'; 
        if find(FG, "CG") eq 0 then delete;
        num_gene = count(primary_sex_det_genes, '|') + 1;
        num = index(FG, '_');
        if num >0 then cgnumber = substr(FG, 1, num-1);
        else cgnumber = FG;
        keep FG primary_sex_det_genes num_gene cgnumber;
        run;

/* Merge on gene symbol information to get some biological understanding */
    proc sort data=neigh2;
        by cgnumber;
        run;

    proc sort data=DMEL548.symbol2proteincg nodupkey out=melclean;
        by cgnumber;
        run;

    data merged;
        merge melclean (in=in1) neigh2 (in=in2) ;
        by cgnumber;
        if in2;
        keep FG symbol primary_sex_det_genes num_gene;
        run;

    proc sort data=merged;
        by DESCENDING num_gene;
        run;

/* Create permanant dataset */
data SEM.dsrp_ggn_primary_neigh_gene_list;
    set merged;
    run;

proc export data=SEM.dsrp_ggn_primary_neigh_gene_list outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_primary_neighborhood_gene_list.csv' dbms=csv replace;
putnames=yes;
run;
