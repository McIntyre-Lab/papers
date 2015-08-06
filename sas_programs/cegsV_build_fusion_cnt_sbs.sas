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
    by Fbgn fusion_id;
    run;

data fusions2;
    set fusions;
    count + 1;
    by Fbgn fusion_id;
    if first.Fbgn then count = 1;
    run;

proc sort data=fusions2;
    by fusion_id;
    run;

/* Pull in sex det corrected symbols */
    data sex;
        set SEM.cegsV_sex_det_stack;
        keep fusion_id symbol_cat;
        run;

    proc sort data=sex nodupkey;
        by symbol_cat fusion_id;
        run;

    data sex2;
        set sex;
        scount + 1;
        by symbol_cat fusion_id;
        if first.symbol_cat then scount = 1;
        run;

    proc sort data=sex2 nodupkey;
        by fusion_id;
        run;

/* For sex Det gene replace fbgn with corrected gene symbol then combine with fusion_id */
    data merged;
        merge fusions2 (in=in1) sex2 (in=in2);
        by fusion_id;
        if in2 then symbol = trim(symbol_cat) || '_' || strip(scount);
        else symbol = trim(FBgn) || '_' || strip(count);
        drop FBgn symbol_cat count scount;
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

    data SEM.cegsV_by_fusion_cnt_sbs;
        set flip;
        drop _name_;
        run;

    data SEM.cegsV_by_fusion_sex_det_cnt_sbs;
        set flip;
        drop _name_ FBgn:;
        run;

/* Clean up */
proc datasets nolist;
    delete flip;
    delete fusions;
    delete fusions2;
    delete merged;
    delete merged2;
    delete sasmacr;
    delete sex;
    delete sex2;
    run; quit;
