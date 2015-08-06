/********************************************************************************
* Create and export dataset for gaussian graphical models using the sex
* determination subset. Then run GGM R script. GGM allows you to select
* edges by an FDR cutoff or by specifiying a number of edges. Two
* different FDR cutoffs were used (0.2, 0.8) along with outputing the
* top 20 edges. GGM was run on the entire DSRP dataset.
********************************************************************************/

data w_sample;
    retain sample;
    format sample $ 20.;
    set SEM.dsrp_sbs_gene_level_sym;
    sample = strip(patRIL) || '_' || strip(matRIL);
    drop patRIL matRIL;
    run;

* I could not get this to export all of the columns (7424). So I ended up
* exporting the file to a csv from JMP genomics. When saving under options select csv.
* saved as /home/jfear/tmp/dsrp_by_gene_sbs.txt;


/* Run GGM On Mated and Virgin */
data _null_;
    call system('Rscript $MCLAB/cegs_sem_sd_paper/r_programs/dsrp_ggm_gene.R');
    run;


