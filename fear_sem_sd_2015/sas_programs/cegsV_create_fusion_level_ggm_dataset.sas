/********************************************************************************
* Export dataset for gaussian graphical models using the sex determination
* subset. Then run GGM R script. GGM allows you to select edges by an FDR
* cutoff or by specifiying a number of edges. Two different FDR cutoffs were
* used (0.2, 0.8) along with outputing the top 20 edges. GGM was run on the
* entire DSRP dataset.
********************************************************************************/
* I could not get this to export all of the columns (32723) in SEM.cegsV_by_fusion_sbs. So I ended up
* exporting the file to a csv from JMP genomics. When saving under options select csv.
*
* Move file to /tmp/data_for_ggm.csv;

/* Run GGM On Mated and Virgin */
data _null_;
    call system('Rscript $MCLAB/cegs_sem_sd_paper/r_programs/cegsV_ggm_fusions.R');
    run;
