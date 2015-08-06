/********************************************************************************
* Create a subset of only the sex determination genes.
********************************************************************************/

/* Pull out Sex Det genes */
    * Using a where statment so that I can keep multi gene fusions as well.;
    data sex;
        set SEM.cegs_virgin_norm_cent;
        where 
        symbol_cat ? 'B52' or 
        symbol_cat ? 'dsx' or 
        symbol_cat ? 'fl(2)d' or 
        symbol_cat ? 'fru' or 
        symbol_cat ? 'dsx' or 
        symbol_cat ? 'her' or 
        symbol_cat ? 'ix' or 
        symbol_cat ? 'msl_2' or 
        symbol_cat ? 'mub' or 
        symbol_cat ? 'ps'  or
        symbol_cat ? 'Psi' or 
        symbol_cat ? 'Rbp1' or 
        symbol_cat ? 'Rm62' or 
        symbol_cat ? 'snf' or 
        symbol_cat ? 'Spf45' or 
        symbol_cat ? 'sqd' or 
        symbol_cat ? 'Sxl' or 
        symbol_cat ? 'tra' or 
        symbol_cat ? 'tra2' or 
        symbol_cat ? 'vir' or 
        symbol_cat ? 'Yp1' or 
        symbol_cat ? 'Yp2' or 
        symbol_cat ? 'Yp3'
        ;
        run;

/* Remove genes not in sex det */
    * Because I used the where statement, there are some off taget genes. I
    * want to remove those genes that have nothing to do with sex det.;
    proc freq data=sex;
        table symbol_cat;
        run; 
        
    data sex2;
        set sex;
        if symbol_cat eq '14-3-3epsilon' then delete;
        if symbol_cat eq 'Aps' then delete;
        if symbol_cat eq 'CG14057|Papst2' then delete;
        if symbol_cat eq 'CG3652|Tps1' then delete;
        if symbol_cat eq 'CG4329|Vps20' then delete;
        if symbol_cat eq 'CG4553|Vps33B' then delete;
        if symbol_cat eq 'CG6044|ventrally-expressed-protein-D' then delete;
        if symbol_cat eq 'CG7741|Fpps' then delete;
        if symbol_cat eq 'CG8129|Fps85D' then delete;
        if symbol_cat eq 'CG9548|epsilonCOP' then delete;
        if symbol_cat eq 'Caps' then delete;
        if symbol_cat eq 'Cpsf100|Slu7' then delete;
        if symbol_cat eq 'Cpsf160' then delete;
        if symbol_cat eq 'Cpsf73' then delete;
        if symbol_cat eq 'Cpsf73|CstF-64' then delete;
        if symbol_cat eq 'Cyp6t2Psi' then delete;
        if symbol_cat eq 'Cyp9f3Psi' then delete;
        if symbol_cat eq 'DNApol-epsilon' then delete;
        if symbol_cat eq 'Dgkepsilon' then delete;
        if symbol_cat eq 'Dgkepsilon|Nacalpha' then delete;
        if symbol_cat eq 'Eps-15' then delete;
        if symbol_cat eq 'Fpps' then delete;
        if symbol_cat eq 'Fps85D' then delete;
        if symbol_cat eq 'Gycbeta100B|stops' then delete;
        if symbol_cat eq 'Lapsyn' then delete;
        if symbol_cat eq 'MP1|vps24' then delete;
        if symbol_cat eq 'Nipsnap' then delete;
        if symbol_cat eq 'Optix' then delete;
        if symbol_cat eq 'Papss' then delete;
        if symbol_cat eq 'Papst2' then delete;
        if symbol_cat eq 'Rbp1-like' then delete;
        if symbol_cat eq 'Reps' then delete;
        if symbol_cat eq 'Rm62' then delete;
        if symbol_cat eq 'Scgbeta|pps' then delete;
        if symbol_cat eq 'Spps|cav' then delete;
        if symbol_cat eq 'Sps2' then delete;
        if symbol_cat eq 'Swim|Vps35' then delete;
        if symbol_cat eq 'Tps1' then delete;
        if symbol_cat eq 'Vps13' then delete;
        if symbol_cat eq 'Vps13|blow' then delete;
        if symbol_cat eq 'Vps16A' then delete;
        if symbol_cat eq 'Vps16B|ca' then delete;
        if symbol_cat eq 'Vps25|beta3GalTII' then delete;
        if symbol_cat eq 'Vps26' then delete;
        if symbol_cat eq 'Vps28' then delete;
        if symbol_cat eq 'Vps33B' then delete;
        if symbol_cat eq 'Vps35' then delete;
        if symbol_cat eq 'Vps36' then delete;
        if symbol_cat eq 'Vps4' then delete;
        if symbol_cat eq 'Vps45' then delete;
        if symbol_cat eq 'cal1|cher' then delete;
        if symbol_cat eq 'calypso|clu' then delete;
        if symbol_cat eq 'caps' then delete;
        if symbol_cat eq 'cher' then delete;
        if symbol_cat eq 'eIF2B-epsilon' then delete;
        if symbol_cat eq 'eIF2B-epsilon|mip130' then delete;
        if symbol_cat eq 'heix' then delete;
        if symbol_cat eq 'l(3)psg2' then delete;
        if symbol_cat eq 'l(3)psg2|spo' then delete;
        if symbol_cat eq 'msps' then delete;
        if symbol_cat eq 'pix' then delete;
        if symbol_cat eq 'pps' then delete;
        if symbol_cat eq 'psh' then delete;
        if symbol_cat eq 'psidin' then delete;
        if symbol_cat eq 'psq' then delete;
        if symbol_cat eq 'pst' then delete;
        if symbol_cat eq 'qkr58E-3|ventrally-expressed-protein-D' then delete;
        if symbol_cat eq 'raps' then delete;
        if symbol_cat eq 'snoRNA:Psi18S-531' then delete;
        if symbol_cat eq 'snoRNA:Psi28S-291' then delete;
        if symbol_cat eq 'stops' then delete;
        if symbol_cat eq 'stumps' then delete;
        if symbol_cat eq 'tral' then delete;
        if symbol_cat eq 'vir-1' then delete;
        if symbol_cat eq 'vps2' then delete;
        if symbol_cat eq 'vps24' then delete;
        if symbol_cat eq 'yps' then delete;
        run;

/* Rename multi-gene fusions */
    * A handful of the fusions pulled have multi gene names. I will rename
    * these so that they only correspond to sex det genes.;
    proc freq data=sex2;
        table symbol_cat;
        run;

    data sex3;
        set sex2;
        if symbol_cat eq 'CG10841|sqd' then symbol_cat = 'sqd';
        if symbol_cat eq 'CG34424|vir' then symbol_cat = 'vir';
        if symbol_cat eq 'Dek|Psi' then symbol_cat = 'Psi';
        if symbol_cat eq 'blos1|tra2' then symbol_cat = 'tra2';
        if symbol_cat eq 'dsx|lds' then symbol_cat = 'dsx';
        if symbol_cat eq 'l(3)73Ah|tra' then symbol_cat = 'tra';
        if symbol_cat eq 'fl(2)d' then symbol_cat = 'fl_2_d';
        drop genes_per_fusion;
        run;

/* Create Permanent sex det stacked dataset */
    data SEM.cegsV_sex_det_stack;
        set sex3;
        run;

    * Check: Look and see which genes are in the dataset;
    /* 
        proc sort data=SEM.cegsV_sex_det_stack;
            by symbol_cat;
            run;
        
        proc freq data=SEM.cegsV_sex_det_stack;
            table symbol_cat;
            run;

        * Genes Present: B52 Psi Rbp1 Spf45 Sxl Yp1 Yp2 Yp3 dsx fl_2_d fru her ix mub ps snf sqd tra tra2 vir;

    */

/* Create a Permanent Sex Det Flag */
    data sex4;
        set sex3;
        keep fusion_id;
        run;

    proc sort data=sex4;
        by fusion_id;
        run;

    proc sort data=DMEL551.fb551_si_fusions_unique;
        by fusion_id;
        run;

    data SEM.cegs_flag_sex_det;
        merge sex4 (in=in1) DMEL551.fb551_si_fusions_unique (in=in2);
        by fusion_id;
        if in1 then flag_sex_det = 1;
        else flag_sex_det = 0;
        keep fusion_id flag_sex_det;
        run;
        
/* Clean-up */
proc datasets nolist;
    delete sex;
    delete sex2;
    delete sex3;
    delete sex4;
    delete merged;
    delete merged2;
    delete fusions;
    run;quit;
