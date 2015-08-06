/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create smaple name based on RIL combination */
    data dsrp;
        retain sample;
        set SEM.dsrp_sbs_combine_sym;
        sample = triM(strip(patRIL)) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* Transpose dataset, because MMC expects genes as rows and obs as cols */
    proc transpose data=dsrp out=flip;
        var :;
        id sample;
        run;

/* Export dataset */
    proc export data=flip outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/mmc/dsrp_data_for_mmc.csv' dbms=csv replace;
        putnames=yes;
        run;

/* output only sex det genes */
    data sex;
        set flip;
        if _name_ eq 'Spf45' or
           _name_ eq 'Sxl_A' or
           _name_ eq 'Sxl_B' or
           _name_ eq 'Sxl_C' or
           _name_ eq 'Yp1' or
           _name_ eq 'Yp2' or
           _name_ eq 'Yp3' or
           _name_ eq 'fl_2_d_A' or
           _name_ eq 'fl_2_d_B' or
           _name_ eq 'fru' or
           _name_ eq 'her' or
           _name_ eq 'ix' or
           _name_ eq 'snf' or
           _name_ eq 'tra' or
           _name_ eq 'tra2_A' or
           _name_ eq 'tra2_B' or
           _name_ eq 'vir';
       run;

    proc export data=sex outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/mmc/dsrp_isoform_sex_det_for_mmc.csv' dbms=csv replace;
        putnames=yes;
        run;
