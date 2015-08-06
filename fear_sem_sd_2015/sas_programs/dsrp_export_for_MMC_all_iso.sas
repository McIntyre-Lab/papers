/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create smaple name based on RIL combination */
    data dsrp;
        retain sample;
        set SEM.dsrp_sbs_symbol;
        sample = triM(strip(patRIL)) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* Transpose dataset, because MMC expects genes as rows and obs as cols */
    proc transpose data=dsrp out=flip;
        var :;
        id sample;
        run;

/* output only sex det genes */
    data sex;
        set flip;
        if prxmatch('/^CG/',_name_) then delete;
       run;

    proc export data=sex outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/mmc/dsrp_all_isoform_sex_det_for_mmc.csv' dbms=csv replace;
        putnames=yes;
        run;
