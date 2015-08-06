/********************************************************************************
* Multicollinearity is going to be a large problem. It will be best to
* collapse isoforms if they have a large correlation coefficient (>=0.8).
* 
* This script calculates the correlation between isoforms and flags them if
* they have a correlation >=0.8.
********************************************************************************/

/* Use a python script to combine highly correlated isoforms within gene */
    %macro combine_by_corr(isoform);
        %if &isoform eq ps %then %do;
            data subset;
            set SEM.dsrp_sex_det_sbs_symbol;
            keep patRIL matRIL &isoform.;
            run;
        %end;
        %else %if &isoform eq tra %then %do;
            data subset;
            set SEM.dsrp_sex_det_sbs_symbol;
            keep patRIL matRIL &isoform._:;
            run;
        %end;
        %else %do;
            data subset;
            set SEM.dsrp_sex_det_sbs_symbol;
            keep patRIL matRIL &isoform.:;
            run;
        %end;

        proc export data=subset outfile='/tmp/corr_subset.csv' dbms=csv replace;
        putnames=yes;
        run;

        data _null_;
            command = 'python $MCLAB/cegs_sem_sd_paper/scripts/clusterByCorr.py -c 0.75 -i /tmp/corr_subset.csv -o $MCLAB/cegs_sem_sd_paper/analysis_output/correlation/' || "&isoform" ||'.csv -t $MCLAB/cegs_sem_sd_paper/analysis_output/correlation/' || "&isoform" ||'_iso_table.csv -g $MCLAB/cegs_sem_sd_paper/analysis_output/correlation/' || "&isoform" ||'.log';
            call system(command);
            run;
    %mend;

    %combine_by_corr(B52);
    %combine_by_corr(fl_2_d);
    %combine_by_corr(fru);
    %combine_by_corr(her);
    %combine_by_corr(ix);
    %combine_by_corr(mub);
    %combine_by_corr(ps);
    %combine_by_corr(Psi);
    %combine_by_corr(Rbp1);
    %combine_by_corr(Rm62);
    %combine_by_corr(snf);
    %combine_by_corr(Spf45);
    %combine_by_corr(sqd);
    %combine_by_corr(Sxl);
    %combine_by_corr(tra);
    %combine_by_corr(tra2);
    %combine_by_corr(vir);
    %combine_by_corr(Yp1);
    %combine_by_corr(Yp2);
    %combine_by_corr(Yp3);


/* Merge all of the genes back together into a single dataset */
    data combined;
        set SEM.dsrp_sex_det_sbs_symbol;
        keep patRIL matRIL;
        run;

    proc sort data=combined;
        by patRIL matRIL;
        run;

    %macro import_merge(isoform);
        proc import datafile="!MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform..csv" out=newIso dbms=csv replace;
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
    %import_merge(B52);
    %import_merge(fl_2_d);
    %import_merge(fru);
    %import_merge(her);
    %import_merge(ix);
    %import_merge(mub);
    %import_merge(ps);
    %import_merge(Psi);
    %import_merge(Rbp1);
    %import_merge(Rm62);
    %import_merge(snf);
    %import_merge(Spf45);
    %import_merge(sqd);
    %import_merge(Sxl);
    %import_merge(tra);
    %import_merge(tra2);
    %import_merge(vir);
    %import_merge(Yp1);
    %import_merge(Yp2);
    %import_merge(Yp3);

    data SEM.dsrp_sex_det_sbs_combine_sym;
        set combined;
        run;


/* Clean Up */
    proc datasets nolist;
        delete subset ;
        delete merged ;
        delete combined ;
        delete newISO;
        run; quit;
