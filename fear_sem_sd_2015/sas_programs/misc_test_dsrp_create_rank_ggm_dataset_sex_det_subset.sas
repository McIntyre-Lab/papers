/********************************************************************************
* Create and export dataset for gaussian graphical models using the sex
* determination subset. Then run GGM R script. GGM allows you to select
* edges by an FDR cutoff or by specifiying a number of edges. Two
* different FDR cutoffs were used (0.2, 0.8) along with outputing the
* top 20 edges. Also Isoforms and gene level variables are used.
********************************************************************************/

/* GGM on Genes in the Sex Det Pathway isoform level */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sexdet_comb_sym_from_rank;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

    proc export data=w_sample outfile='/tmp/data_for_ggm.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Run R-script;
    data _null_;
        call system('Rscript /home/jfear/mclab/cegs_sem_sd_paper/r_programs/dsrp_rank_ggm_sex_det_subset.R');
        run;

/* Clean up */
proc datasets nolist;
    delete w_sample;
    run; quit;
