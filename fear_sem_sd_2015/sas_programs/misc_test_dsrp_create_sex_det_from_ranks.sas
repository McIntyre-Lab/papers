/* Create Initial Sex Det dataset from Ranks*/
    data sex;
        set SEM.dsrp_stack_rank;
        if 
        symbol eq 'B52' or 
        symbol eq 'dsx' or 
        symbol eq 'fl(2)d' or 
        symbol eq 'fru' or 
        symbol eq 'her' or 
        symbol eq 'ix' or 
        symbol eq 'msl-2' or 
        symbol eq 'mub' or 
        symbol eq 'ps'  or
        symbol eq 'Psi' or 
        symbol eq 'Rbp1' or 
        symbol eq 'Rm62' or 
        symbol eq 'snf' or 
        symbol eq 'Spf45' or 
        symbol eq 'sqd' or 
        symbol eq 'Sxl' or 
        symbol eq 'tra' or 
        symbol eq 'tra2' or 
        symbol eq 'vir' or 
        symbol eq 'Yp1' or 
        symbol eq 'Yp2' or 
        symbol eq 'Yp3'
        ;
        run;

    proc freq data=sex;
        table symbol;
        run;

    * Create Permanent sex det stacked dataset;
    data SEM.dsrp_sexdet_stack_from_rank;
        set sex;
        run;

    * Create Permanent Sex Det side-by-side with symbol;
    proc sort data=sex;
        by patRIL matRIL symbol cgnumber;
        run;

    * Append a counter to symbol so I can keep isofroms separate;
    data counts;
        set sex;
        count + 1;
        by patRIL matRIL symbol;
        if first.symbol then count = 1;
        run;

    data counts2;
        set counts;
        sym = trim(symbol) || '_' || strip(count);
        run;

    proc sort data=counts2;
        by patRIL matRIL;
        run;

    proc transpose data=counts2 out=cgflip2;
        by patRIL matRIL;
        var expRank;
        id sym;
        run;

    data SEM.dsrp_sexdet_sbs_sym_from_rank;
        set cgflip2;
        drop _name_ _label_;
        run;

    * Clean-up;
    proc datasets nolist;
        delete sex;
        delete counts;
        delete counts2;
        delete cgflip2;
        run;quit;

/* Calculate Correlation to identify highly correlated (>=0.8) isoforms and average ranks */
    %macro make_corr(isoform);
        proc corr data=SEM.dsrp_sexdet_sbs_sym_from_rank out=corr noprint;
            var &isoform:;
            run;

        data corr2;
            set corr;
            where _type_ eq 'CORR';
            rename _name_ = iso1;
            label _name_ = ' ';
            run;

        proc sort data=corr2;
            by iso1;
            run;

        proc transpose data=corr2 out=flip;
            by iso1;
            var &isoform:;
            run;
            
        data flip2;
            set flip;
            rename _name_ = iso2;
            label _name_ = ' ';
            if col1 > 0.8 then flag_corr = 1; else flag_corr = 0;
            drop col1;
            run;

        proc transpose data=flip2 out=flip3;
            by iso1;
            var flag_corr;
            id iso2;
            run;
            
        proc sort data=flip3 out=flipsort;
            by &isoform._1;
            run;

        data &isoform;
            set flipsort;
            drop _name_;
            run;

        proc datasets nolist;
            delete corr;
            delete corr2;
            delete flip;
            delete flip2;
            delete flip3;
            delete flipsort;
            run; quit;

        options ls=200;
        proc print data=&isoform;
        run;

        proc datasets nolist;
            delete &isoform;
            run; quit;

    %mend;

    %make_corr(B52);
    %make_corr(fl_2_d);
    %make_corr(fru);
    %make_corr(her);
    %make_corr(ix);
    %make_corr(mub);
    %make_corr(ps);
    %make_corr(Psi);
    %make_corr(Rbp1);
    %make_corr(Rm62);
    %make_corr(snf);
    %make_corr(Spf45);
    %make_corr(sqd);
    %make_corr(Sxl);
    %make_corr(tra);
    %make_corr(tra2);
    %make_corr(vir);
    %make_corr(Yp1);
    %make_corr(Yp2);
    %make_corr(Yp3);

    data SEM.dsrp_sexdet_comb_sym_from_rank;
        set SEM.dsrp_sexdet_sbs_sym_from_rank;

        B52_A = mean(B52_4,B52_5,B52_6,B52_7);
        B52_B = mean(B52_1,B52_2,B52_3,B52_8,B52_9,B52_10);

        Psi = Psi_1;

        Rdp1 = mean(Rbp1_1,Rbp1_2);

        Rm62 = Rm62_1;

        Spf45_A = Spf45_1;
        Spf45_B = mean(Spf45_2,Spf45_3);

        Sxl_A = mean(Sxl_1,Sxl_3,Sxl_5,Sxl_6,Sxl_7,Sxl_9,Sxl_10,Sxl_11,Sxl_13,Sxl_15,Sxl_17,Sxl_19,Sxl_20,Sxl_21,Sxl_16,Sxl_22,Sxl_23);
        Sxl_B = mean(Sxl_2,Sxl_4,Sxl_8,Sxl_14);
        Sxl_C = mean(Sxl_12,Sxl_18,Sxl_24);

        Yp1 = Yp1_1;

        Yp2 = Yp2_1;

        Yp3 = Yp3_1;

        fl_2_d_A = fl_2_d_1;
        fl_2_d_B = mean(fl_2_d_2,fl_2_d_3,fl_2_d_4);

        fru = mean(fru_1,fru_2,fru_3,fru_4,fru_5,fru_6,fru_7,fru_8,fru_9,fru_10,fru_11,fru_12,fru_13,fru_14,fru_15);

        her = her_1;

        ix = ix_1;

        mub_A = mean(mub_3,mub_9);
        mub_B = mean(mub_4,mub_5,mub_8);
        mub_c = mean(mub_1,mub_2,mub_6,mub_7,mub_10,mub_11);

        ps = ps_1;

        snf = snf_1;

        sqd_A = mean(sqd_1,sqd_2,sqd_3,sqd_5);
        sqd_B = sqd_4;

        tra = tra_1;

        tra2_A = mean(tra2_1,tra2_3,tra2_4,tra2_5,tra2_6,tra2_7);
        tra2_B = tra2_2;

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

