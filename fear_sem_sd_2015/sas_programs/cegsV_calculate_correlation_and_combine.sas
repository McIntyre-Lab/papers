/********************************************************************************
* Multicollinearity is going to be a large problem. It will be best to
* collapse isoforms if they have a large correlation coefficient (>=0.8).
* 
* This script calculates the correlation between isoforms and flags them if
* they have a correlation >=0.8.
********************************************************************************/
/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

* Copy locally for faster processing;
    data cegsV_sbs_symbol;
        set SEM.cegsV_by_fusion_cnt_sbs;
        run;

    data gene_list;
        set SEM.cegsV_gene_list;
        run;

/* Use a python script to combine highly correlated fusions within gene */
    %macro combine_by_corr(gene);
        %if %varexist(cegsV_sbs_symbol, &gene) %then %do;
            data subset;
                set cegsV_sbs_symbol;
                keep line &gene.;
                run;
        %end;
        %else %do;
            data subset;
                set cegsV_sbs_symbol;
                keep line &gene._:;
                run;
        %end;

        proc export data=subset outfile='/tmp/corr_subset.csv' dbms=csv replace;
            putnames=yes;
            run;

        data _null_;
            command = 'python $MCLAB/cegs_sem_sd_paper/scripts/cegs_clusterByCorr.py -c 0.75 -i /tmp/corr_subset.csv -o $HOME/tmp/analysis_output/correlation/' || "&gene" ||'.csv -t $HOME/tmp/analysis_output/correlation/' || "&gene" ||'_fusion_table.csv';
            call system(command);
            run;
    %mend;

    proc printto log='$HOME/tmp/combine.log' new; run;
    %iterdataset(dataset=gene_list, function=%nrstr(%combine_by_corr(&symbol);));

/* Merge all of the genes back together into a single dataset */
    data combined;
        set cegsV_sbs_symbol;
        keep line;
        run;

    proc sort data=combined;
        by line;
        run;

    %macro import_merge(gene);
        proc import datafile="!HOME/tmp/analysis_output/correlation/&gene..csv" out=newIso dbms=csv replace;
            getnames=yes;
            run;

        proc sort data=newISO;
            by line;
            run;

        data merged;
            merge combined newISO;
            by line;
            run;

        data combined;
            set merged;
            run;
    %mend;

    %iterdataset(dataset=gene_list, function=%nrstr(%import_merge(&symbol);));
    proc printto log=log new; run;

    data SEM.cegsV_by_fusion_comb_sbs;
        set combined;
        run;

/* Clean Up */
    proc datasets nolist;
        delete subset ;
        delete merged ;
        delete combined ;
        delete newISO;
        run; quit;
