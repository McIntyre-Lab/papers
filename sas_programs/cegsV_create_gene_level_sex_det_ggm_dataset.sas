/********************************************************************************
* Export dataset for gaussian graphical models using the sex determination
* subset. Then run GGM R script. GGM allows you to select edges by an FDR
* cutoff or by specifiying a number of edges. Two different FDR cutoffs were
* used (0.2, 0.8) along with outputing the top 20 edges. GGM was run on the
* entire DSRP dataset.
********************************************************************************/
* I could not get this to export all of the columns (8823) in SEM.cegsV_by_gene_sbs. So I ended up
* exporting the file to a csv from JMP genomics. When saving under options select csv.
*
* Saved SEM.cegsV_by_gene_sbs /home/jfear/tmp/cegsv_by_gene_sbs.txt;

proc export data=SEM.cegsV_by_gene_sex_det_sbs outfile='/home/jfear/tmp/cegsv_by_gene_sex_det_sbs.txt' dbms=csv replace;
    putnames =yes;
    run;

/* Run GGM on CEGS Virgin */
data _null_;
    call system('Rscript $MCLAB/cegs_sem_sd_paper/r_programs/cegsV_ggm_gene_sex_det.R');
    run;
