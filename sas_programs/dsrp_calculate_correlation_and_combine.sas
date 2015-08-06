/********************************************************************************
* Multicollinearity is going to be a large problem. It will be best to
* collapse isoforms if they have a large correlation coefficient (>=0.8).
* 
* This script calculates the correlation between isoforms and flags them if
* they have a correlation >=0.8.
********************************************************************************/

libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Create Gene list */
data dsrp_stack_short;
    set SEM.dsrp_stack;
    keep patRIL matRIL cgnumber;
    run;

data gene_list;
    set dsrp_stack_short;
    count = index(cgnumber, '_');
    if count >0 then do;
        symbol = trim(substr(cgnumber,1, count-1));
    end;
    else symbol = cgnumber;
    keep symbol;
    run;

proc sort data=gene_list nodupkey;
    by symbol;
    run; * 7422 obs;

/* Use a python script to combine highly correlated isoforms within gene */
    %macro combine_by_corr(isoform);
        data subset;
            set SEM.dsrp;
            keep patRIL matRIL &isoform.:;
            run;

        proc export data=subset outfile='/tmp/corr_subset.csv' dbms=csv replace;
            putnames=yes;
            run;

        data _null_;
            command = 'python $MCLAB/cegs_sem_sd_paper/scripts/clusterByCorr.py -c 0.75 -i /tmp/corr_subset.csv -o $HOME/tmp/analysis_output/correlation/' || "&isoform" ||'_dsx.csv -t $HOME/tmp/analysis_output/correlation/' || "&isoform" ||'_iso_table_dsx.csv';
            call system(command);
            run;
    %mend;

    %iterdataset(dataset=gene_list, function=%nrstr(%combine_by_corr(&symbol);));


/* Merge all of the genes back together into a single dataset */
    data combined;
        set SEM.dsrp;
        keep patRIL matRIL;
        run;

    proc sort data=combined;
        by patRIL matRIL;
        run;

    %macro import_merge(isoform);
        proc import datafile="!HOME/tmp/analysis_output/correlation/&isoform._dsx.csv" out=newIso dbms=csv replace;
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

    data SEM.dsrp_combine;
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
