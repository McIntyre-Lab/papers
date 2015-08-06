/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create local copy of main dataset and split FBgn_cat up */
data fusions;
    set SEM.cegs_virgin_norm_cent;
    num = index(Fbgn_cat, '|');
    if num > 0 then Fbgn = substr(Fbgn_cat, 1, num - 1);
    else Fbgn = Fbgn_cat;
    keep fusion_id FBgn;
    run;

proc sort data=fusions nodupkey;
    by fusion_id;
    run;

/* Pull in sex det corrected symbols */
    data sex;
        set SEM.cegsV_sex_det_stack;
        keep fusion_id symbol_cat;
        run;

    proc sort data=sex nodupkey;
        by fusion_id;
        run;

/* For sex Det gene replace fbgn with corrected gene symbol then combine with fusion_id */
    data merged;
        merge fusions (in=in1) sex (in=in2);
        by fusion_id;
        if in2 then symbol = trim(symbol_cat) || '_' || trim(fusion_id);
        else symbol = trim(FBgn) || '_' || trim(fusion_id);
        drop FBgn symbol_cat;
        run;

/* Merge to the main dataset */
    proc sort data=merged;
        by fusion_id;
        run;

    proc sort data=SEM.cegs_virgin_norm_cent;
        by fusion_id;
        run;

    data merged2;
        merge merged (in=in1) SEM.cegs_virgin_norm_cent;
        by fusion_id;
        keep symbol line mean_exp;
        run;

/* Flip and create permanent dataset */
    proc sort data=merged2;
        by line;
        run;

    proc transpose data=merged2 out=flip;
        by line;
        var mean_exp;
        id symbol;
        run;

    data SEM.cegsV_by_fusion_sbs;
        set flip;
        drop _name_;
        run;

    data SEM.cegsV_by_fusion_sex_det_sbs;
        set flip;
        drop _name_ FBgn:;
        run;
