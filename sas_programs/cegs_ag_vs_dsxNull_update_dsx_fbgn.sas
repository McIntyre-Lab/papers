/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Check if there are FBgns in dsxNull that are not in the FB5.51 primary fbgn */
    * copy needed datasets and rename vars so I can merge;
        data fb551;
            set DMEL551.symbol2fbgn;
            run;

        data dsxNull_repressed;
            set SEM.DSXNullF_repressed;
            rename FBgn_cat = primary_fbgn;
            run;

        data dsxNull_induced;
            set SEM.DSXNullF_induced;
            rename FBgn_cat = primary_fbgn;
            run;

    * Merge data and see if there are FBGNs that don't match;
        proc sort data=fb551;
            by primary_fbgn;
            run;

        proc sort data=dsxNull_repressed;
            by primary_fbgn;
            run;

        proc sort data=dsxNull_induced;
            by primary_fbgn;
            run;

        data merged oops;
            merge fb551 (in=in1) dsxNull_repressed (in=in2) dsxNull_induced (in=in3);
            by primary_fbgn;
            if in2 then flag_repressed = 1;
            else flag_repressed = 0;
            if in3 then flag_induced = 1;
            else flag_induced = 0;
            if in1 and (in2 or in3) then output merged;
            else if in2 or in3 then output oops;
            run; * 23 obs in oops :-( ;

/* Import FBgn conversions to see if I can resolve these 16 */
     data work.fbgn2;
         infile '!MCLAB/useful_dmel_data/flybase551/flybase_files/fbgn_annotation_ID_fb_2012_06.tsv' delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=6;
            informat symbol $81. ;
            informat primary_fbgn $13. ;
            informat secondary_fbgn $131. ;
            informat primary_cg $13. ;
            informat secondary_cg $26. ;
            format symbol $81. ;
            format primary_fbgn $13. ;
            format secondary_fbgn $131. ;
            format primary_cg $13. ;
            format secondary_cg $26. ;
         input
                     symbol $
                     primary_fbgn $
                     secondary_fbgn $
                     primary_cg $
                     secondary_cg $
         ;
         run;

    data fbgn_all; 
        retain primary_fbgn;
        length fbgn $160.;
        set fbgn2;
        if secondary_fbgn = ' ' then fbgn = primary_fbgn;
        else fbgn = trim(primary_fbgn) || ',' || trim(secondary_fbgn);
        keep primary_fbgn fbgn;
        run;


/* Iterate through DSX genelist and add corresponding Fbgn form Fb5.51 */
    data dsxFb551_repressed;
        set dsxnull_repressed;
        rename primary_fbgn = fb530;
        run;

    data dsxFb551_induced;
        set dsxnull_induced;
        rename primary_fbgn = fb530;
        run;

    %macro iter_fbgn(tab, fbgn);
        data f551;
            set fbgn_all;
            where fbgn ? "&fbgn";
            run;

        data _null_;
            set f551;
            call symput("primary", primary_fbgn);
            run;

        data &tab;
            set &tab;
            if "&fbgn" = "FBgn0261701" then do;
                if fb530 = "&fbgn" then fb551 = "FBgn0262838";
            end;
            else do;
                if fb530 = "&fbgn" then fb551 = "&primary";
            end;
            run;
    %mend;

    %iterdataset(dataset=dsxNull_repressed,function=%nrstr(%iter_fbgn(dsxFb551_repressed, &primary_fbgn);));
    %iterdataset(dataset=dsxNull_induced,function=%nrstr(%iter_fbgn(dsxFb551_induced, &primary_fbgn);));

    * Everything worked perfectly except for FBgn0261701;
    * I pulled the gene symbol from the FB5.30 data using grep "CG42737" which
    * corresponds to "FBgn0262838" in FB5.51;

/* make dataset */
    data SEM.DSXNullF_repressed_w_fb551;
        set dsxFb551_repressed;
        run;

    data SEM.DSXNullF_induced_w_fb551;
        set dsxFb551_induced;
        run;

    * There are two genes in Fb530 that were merged in Fb551 to we only need one copy;
    proc sort data=SEM.DSXNullF_repressed_w_fb551 nodupkey;
        by fb551;
        run;

    proc sort data=SEM.DSXNullF_induced_w_fb551 nodupkey;
        by fb551;
        run;

/* Clean Up */
    proc datasets nolist;
        delete dsxfb551_repressed;
        delete dsxfb551_induced;
        delete dsxnull_repressed;
        delete dsxnull_induced;
        delete f551;
        delete fb551;
        delete fbgn2;
        delete fbgn_all;
        delete merged;
        delete oops;
        run; quit;
