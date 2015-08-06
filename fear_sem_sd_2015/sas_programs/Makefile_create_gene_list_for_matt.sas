libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* DSPR */
    data sex;
        set SEM.dsrp_sex_det_stack;
        num = index(cgnumber, '_');
        if num > 0 then cgnum = substr(cgnumber, 1, num - 1);
        else cgnum = cgnumber;
        keep cgnum symbol FBgn;
        run;

    proc sort data= sex nodupkey;
        by symbol;
        run;

    proc export data=sex outfile='/home/jfear/mclab/cegs_sem_sd_paper/reports/dspr_for_matt.csv' dbms=csv replace;
        putnames=yes;
        run;

    proc datasets nolist;
        delete sex;
        run;

/* Create DSPR RIL list */
    data dspr;
        set SEM.dsrp_sbs_gene_level_sym;
        keep patRIL matRIL;
        run;

    proc sort data= dspr nodupkey;
        by patRIL matRIL;
        run;

    proc export data=dspr outfile='/home/jfear/mclab/cegs_sem_sd_paper/reports/dspr_RIL_ids_for_matt.csv' dbms=csv replace;
        putnames=yes;
        run;



/* CEGS */
    data sex;
        set SEM.cegsV_sex_det_stack;
        if index(FBgn_cat, '|') then delete;
        rename symbol_cat = symbol;
        rename FBgn_cat = Fbgn;
        keep symbol_cat FBgn_cat;
        run;

    proc sort data= sex nodupkey;
        by symbol;
        run;

    proc export data=sex outfile='/home/jfear/mclab/cegs_sem_sd_paper/reports/cegs_for_matt.csv' dbms=csv replace;
        putnames=yes;
        run;

    proc datasets nolist;
        delete sex;
        run;

/* Create CEGS Line List */
    data cegs;
        set SEM.cegs_virgin_norm_cent;
        keep line;
        run;

    proc sort data= cegs nodupkey;
        by line;
        run;

    proc export data=cegs outfile='/home/jfear/mclab/cegs_sem_sd_paper/reports/cegs_line_list_for_matt.csv' dbms=csv replace;
        putnames=yes;
        run;

