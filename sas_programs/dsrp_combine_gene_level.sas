/********************************************************************************
* I will also collapse everything to the gene level. 
********************************************************************************/

/*
proc contents data= SEM.dsrp_sex_det_sbs_symbol varnum;
run;
*/

data SEM.dsrp_sex_det_sbs_gene_level_sym;
    set SEM.dsrp_sex_det_sbs_symbol;

    B52 = mean(B52_4,B52_5,B52_6,B52_7,B52_1,B52_2,B52_3,B52_8,B52_9,B52_10);

    Psi = Psi_1;

    Rdp1 = mean(Rbp1_1,Rbp1_2);

    Rm62 = Rm62_1;

    Spf45 = mean(Spf45_2,Spf45_3,Spf45_1);

    Sxl = mean(Sxl_1,Sxl_3,Sxl_5,Sxl_6,Sxl_7,Sxl_9,Sxl_10,Sxl_11,Sxl_13,Sxl_15,Sxl_17,Sxl_19,Sxl_20,Sxl_21,Sxl_16,Sxl_22,Sxl_23,Sxl_2,Sxl_4,Sxl_8,Sxl_14,Sxl_12,Sxl_18,Sxl_24);

    Yp1 = Yp1_1;

    Yp2 = Yp2_1;

    Yp3 = Yp3_1;

    fl_2_d = mean(fl_2_d_2,fl_2_d_3,fl_2_d_4,fl_2_d_1);

    fru = mean(fru_1,fru_2,fru_3,fru_4,fru_5,fru_6,fru_7,fru_8,fru_9,fru_10,fru_11,fru_12,fru_13,fru_14,fru_15);

    her = her_1;

    ix = ix_1;

    mub = mean(mub_3,mub_9,mub_4,mub_5,mub_8,mub_1,mub_2,mub_6,mub_7,mub_10,mub_11);

    ps = ps_1;

    snf = snf_1;

    sqd = mean(sqd_1,sqd_2,sqd_3,sqd_5,sqd_4);

    tra = tra_1;

    tra2 = mean(tra2_1,tra2_3,tra2_4,tra2_5,tra2_6,tra2_7,tra2_2);

    vir = vir_1;

    drop B52_4 B52_5 B52_6 B52_7 B52_1 B52_2 B52_3 B52_8 B52_9 B52_10 Psi_1
    Rbp1_1 Rbp1_2 Rm62_1 Spf45_1 Spf45_2 Spf45_3 Sxl_1 Sxl_2 Sxl_3 Sxl_4 Sxl_5
    Sxl_6 Sxl_7 Sxl_8 Sxl_9 Sxl_10 Sxl_11 Sxl_12 Sxl_13 Sxl_14 Sxl_15 Sxl_16
    Sxl_17 Sxl_18 Sxl_19 Sxl_20 Sxl_21 Sxl_22 Sxl_23 Sxl_24 Yp1_1 Yp2_1 Yp3_1
    fl_2_d_1 fl_2_d_2 fl_2_d_3 fl_2_d_4 fru_1 fru_2 fru_3 fru_4 fru_5 fru_6
    fru_7 fru_8 fru_9 fru_10 fru_11 fru_12 fru_13 fru_14 fru_15 her_1 ix_1
    mub_1 mub_2 mub_3 mub_4 mub_5 mub_6 mub_7 mub_8 mub_9 mub_10 mub_11 ps_1
    snf_1 sqd_1 sqd_2 sqd_3 sqd_4 sqd_5 tra_1 tra2_1 tra2_2 tra2_3 tra2_4
    tra2_5 tra2_6 tra2_7 vir_1
    ;
    run;

