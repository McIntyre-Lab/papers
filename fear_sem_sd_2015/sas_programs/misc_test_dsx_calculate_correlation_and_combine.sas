/********************************************************************************
* Multicollinearity is going to be a large problem. It will be best to
* collapse isoforms if they have a large correlation coefficient (>=0.8).
********************************************************************************/

/* Create Gene list */
data gene_list;
    set SEM.dsxnullf_repressed_cgnumber;
    count = index(cgnumber, '_');
    symbol = trim(substr(cgnumber,1, count-1));
    keep symbol;
    run;

proc sort data=gene_list nodupkey;
    by symbol;
    run; * 35 obs;

/* Use a python script to combine highly correlated isoforms within gene */
    %macro combine_by_corr(isoform);
        data subset;
            set SEM.dsrp_dsx_sbs_symbol;
            keep patRIL matRIL &isoform.:;
            run;

        proc export data=subset outfile='/tmp/corr_subset.csv' dbms=csv replace;
            putnames=yes;
            run;

        data _null_;
            command = 'python $MCLAB/cegs_sem_sd_paper/scripts/clusterByCorr.py -c 0.75 -i /tmp/corr_subset.csv -o $MCLAB/cegs_sem_sd_paper/analysis_output/correlation/' || "&isoform" ||'_dsx.csv -t $MCLAB/cegs_sem_sd_paper/analysis_output/correlation/' || "&isoform" ||'_iso_table_dsx.csv -g $MCLAB/cegs_sem_sd_paper/analysis_output/correlation/' || "&isoform" ||'_dsx.log';
            call system(command);
            run;
    %mend;

    %iterdataset(dataset=gene_list, function=%nrstr(%combine_by_corr(&symbol);));


/* Merge all of the genes back together into a single dataset */
    data combined;
        set SEM.dsrp_dsx_sbs_symbol;
        keep patRIL matRIL;
        run;

    proc sort data=combined;
        by patRIL matRIL;
        run;

    %macro import_merge(isoform);
        proc import datafile="!MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform._dsx.csv" out=newIso dbms=csv replace;
            getnames=yes;
            run;

        proc sort data=newISO;
            by patRIL matRIL;
            run;

        data merged;
            merge combined newISO;
            by patRIL matRIL;
            run;

        data combined;
            set merged;
            run;
    %mend;
    %iterdataset(dataset=gene_list, function=%nrstr(%import_merge(&symbol);));

    data SEM.dsrp_dsx_sbs_combine_sym;
        set combined;
        run;

/* Clean Up */
    proc datasets nolist;
        delete gene_list ;
        delete subset ;
        delete merged ;
        delete combined ;
        delete newISO;
        run; quit;
