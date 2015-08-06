libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Import Chang DEG */
proc import datafile='!MCLAB/cegs_sem_sd_paper/original_data/Chang2011_TableS10_deg.csv' out=Chang2011 dbms=csv replace;
    getnames=yes;
    run;

data Chang2011;
    set Chang2011;
    if FvM_FDR < 0.05 and FvM_FC <1 then flag_chang_male_bias = 1;
    else flag_chang_male_bias = 0;
    if FvM_FDR < 0.05 and FvM_FC >1 then flag_chang_female_bias = 1;
    else flag_chang_female_bias = 0;

    if FvTra_FDR < 0.2 and FvTra_FC <1 then flag_chang_tra_bias = 1;
    else flag_chang_tra_bias = 0;
    if FvTra_FDR < 0.2 and FvTra_FC >1 then flag_chang_traf_bias = 1;
    else flag_chang_traf_bias = 0;

    if FvM_FDR < 0.05 and FvTra_FDR < 0.2 and ((FvM_FC<1 and FvTra_FC<1) or (FvM_FC>1 and FvTra_FC>1)) then flag_chang_ds_tra = 1;
    else flag_chang_ds_tra = 0;

    if FvM_FDR < 0.05 and FvTra_FDR > 0.2 then flag_chang_ups_tra = 1;
    else flag_chang_ups_tra = 0;

    drop FvM_FDR FvM_FC FvTra_FDR FvTra_FC;
    run;

/* Import Chang tra binding sites */
proc import datafile='!MCLAB/cegs_sem_sd_paper/original_data/Chang2011_TableS20_tra_bs.csv' out=Chang2011_bs dbms=csv replace;
    getnames=yes;
    run;

/* Import Goldman MF deg */
proc import datafile='!MCLAB/cegs_sem_sd_paper/original_data/Goldman2007_TableS1_MF_deg.csv' out=Goldman2007 dbms=csv replace;
    getnames=yes;
    run;

data Goldman2007;
    set Goldman2007;
    if FvsM > 1 then flag_goldman_female_bias = 1;
    else flag_goldman_female_bias = 0;
    if FvsM < 1 then flag_goldman_male_bias = 1;
    else flag_goldman_male_bias = 0;
    drop cgnumber FvsM FvsM_FDR;
    run;

/* Import Goldman DS TRA */
proc import datafile='!MCLAB/cegs_sem_sd_paper/original_data/Goldman2007_TableS2_genes_ds_tra.csv' out=Goldman2007_ds_tra dbms=csv replace;
    getnames=yes;
    run;

data Goldman2007_ds_tra;
    set Goldman2007_ds_tra;
    flag_goldman_ds_tra = 1;
    drop cgnumber FvsM FvsTra;
    run;


/* Merge apaloza */
* Ag merge;
proc sort data=SEM.cegsV_gene2fbgn;
    by gene;
    run;

proc sort data=SEM.cegsV_ag_yp2_added_genes;
    by gene;
    run;

data ag_merge;
    merge SEM.cegsV_gene2Fbgn (in=in1) SEM.cegsV_ag_yp2_added_genes (in=in2);
    by gene;
    if in2 then flag_sem_added_gene = 1;
    else flag_sem_added_gene = 0;
    drop model gene;
    run;

* dsx merge;
data induced;
    set SEM.dsxNullF_induced_w_fb551;
    rename fb551 = primary_fbgn;
    keep fb551;
    run;

data repressed;
    set SEM.dsxNullF_repressed_w_fb551;
    rename fb551 = primary_fbgn;
    keep fb551;
    run;

proc sort data=ag_merge;
    by primary_fbgn;
    run;

proc sort data=induced;
    by primary_fbgn;
    run;

proc sort data=repressed;
    by primary_fbgn;
    run;

data dsx_merge;
    merge ag_merge (in=in1) induced (in=in2) repressed (in=in3);
    by primary_fbgn;
    if in2 then flag_dsxNull_induced = 1;
    else flag_dsxNull_induced = 0;
    if in3 then flag_dsxNull_repressed = 1;
    else flag_dsxNull_repressed = 0;
    run;


* chang merge;
proc sort data=chang2011;
    by primary_fbgn;
    run;

proc sort data=chang2011_bs;
    by primary_fbgn;
    run;

proc sort data=dsx_merge;
    by primary_fbgn;
    run;

data chang_merge;
    merge dsx_merge (in=in1) chang2011 (in=in2) chang2011_bs (in=in3);
    by primary_fbgn;
    if not in1 then delete;
    if not in2 then do;
        flag_chang_male_bias = 0;
        flag_chang_female_bias = 0;
        flag_chang_tra_bias = 0;
        flag_chang_traf_bias = 0;
        flag_chang_ds_tra = 0;
        flag_chang_ups_tra = 0;
    end;
    if in3 then flag_tra_bs = 1;
    else flag_tra_bs = 0;
    drop symbol;
    run;

* goldman Merge ;
proc sort data=Goldman2007;
    by primary_fbgn;
    run;

proc sort data=Goldman2007_ds_tra;
    by primary_fbgn;
    run;

proc sort data=chang_merge;
    by primary_fbgn;
    run;

data goldman_merge;
    merge chang_merge (in=in1) Goldman2007 (in=in2) Goldman2007_ds_tra (in=in3);
    by primary_fbgn;
    if not in1 then delete;
    if not in2 then do;
        flag_goldman_female_bias = 0;
        flag_goldman_male_bias = 0;
    end;
    if not in3 then do;
        flag_goldman_ds_tra = 0;
    end;
    drop symbol;
    run;

* annotation merge;
proc sort data=DMEL551.symbol2fbgn;
    by primary_fbgn;
    run;

proc sort data=goldman_merge;
    by primary_fbgn;
    run;

data for_lauren;
    retain primary_fbgn symbol;
    merge goldman_merge (in=in1) DMEL551.symbol2fbgn (in=in2);
    by primary_fbgn;
    if in1;
    run;

data SEM.validation_set;
    set for_lauren;
    run;

proc contents data=SEM.validation_set;
run;
