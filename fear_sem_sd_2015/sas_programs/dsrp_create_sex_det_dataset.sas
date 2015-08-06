/********************************************************************************
* Merges FlyBase 5.48 gene annotations and creates a subset of the dsrp data
* for just the sex determination genes.
********************************************************************************/

/* Merge on Gene Information */
    * Flip microarray data so I can merge ;
    proc sort data=SEM.dsrp;
        by patRIL matRIL;
        run;

    proc transpose data=SEM.dsrp out=flip;
        by patRIL matRIL;
        var C:;
        run;

    data flip;
        set flip;
        rename _name_ = protein_CG_number;
        label _name_ = ' ';
        run;

    * Merge by protein CG number. Some of the microarray genes have been collapsed
    * to the gene level. So the merge will fail. I will need to merge these by
    * regular CG number.;
    proc sort data=flip;
        by protein_CG_number;
        run;

    proc sort data=DMEL548.symbol2proteincg;
        by protein_CG_number;
        run;

    data mergeOK mergeNO;
        merge flip (in=in1) DMEL548.symbol2proteincg (in=in2);
        by protein_CG_number;
        if in1 and in2 then output mergeOK;
        else if in1 and not in2 then output mergeNO;
        run;

    * Reduce annotations down to gene level for smooth merging;
    data dmel;
        set DMEL548.symbol2proteincg;
        drop FBtr protein_CG_number;
        run;

    proc sort data=dmel nodupkey;
        by cgnumber;
        run;

    * clean up mergeNO dataset prior to merging;
    data mergeNO2;
        set mergeNO;
        drop _LABEL_ symbol FBgn FBtr cgnumber;
        run;

    data mergeNO3;
        set mergeNO2;
        rename protein_CG_number = cgnumber;
        run;

    * merge the mergeNO3 by cgnumber;
    proc sort data=mergeNO3;
        by cgnumber;
        run;

    proc freq data=mergeNO3 noprint;
        table cgnumber /out=freqs;
        run;

    data mergeOK2 oops;
        merge mergeNO3 (in=in1) dmel (in=in2);
        by cgnumber;
        if in1 and in2 then output mergeOK2;
        else if in1 and not in2 then output oops;
        run;

    * Clean mergeNO and Set mergeOks;
    data mergeOKclean;
        set mergeok;
        drop _label_ cgnumber;
        run;

    data mergeOKclean2;
        set mergeOKclean;
        rename protein_CG_number = cgnumber;
        run;

    data merged;
        set mergeOKclean2 mergeOK2;
        run;

/* Save Permanent DSRP stacked dataset */
    data SEM.dsrp_stack;
        set merged;
        run;

/* Pull out Sex Det genes */
    data sex;
        set merged;
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
    data SEM.dsrp_sex_det_stack;
        set sex;
        run;

    * Check: Look and see which genes are in the dataset;
    /* 
    proc sort data=SEM.dsrp_sex_det_stack;
        by symbol;
        run;
    
    proc freq data=SEM.dsrp_sex_det_stack;
        table symbol;
        run;

    * Genes Present: B52 fl(2)d fru her ix mub ps Psi Rbp1 Rm62 snf Spf45 sqd Sxl tra tra2 vir Yp1 Yp2 Yp3;

    */

/* Create Permanent Sex det side-by-side with cgnumber */
    proc sort data=SEM.dsrp_sex_det_stack;
        by patRIL matRIL;
        run;

    proc transpose data=SEM.dsrp_sex_det_stack out=cgflip;
        by patRIL matRIL;
        var col1;
        id cgnumber;
        run;

    data SEM.dsrp_sex_det_sbs_cg;
        set cgflip;
        drop _name_;
        run;

/* Create Permanent Sex Det side-by-side with symbol */
    proc sort data=SEM.dsrp_sex_det_stack;
        by patRIL matRIL symbol cgnumber;
        run;

    * Append a isoform information to symbol so I can keep isofroms separate;
    data suffix;
        set SEM.dsrp_sex_det_stack;
        num = index(cgnumber, '_');
        if num > 0 then do;
            suffix = substr(cgnumber, num+1);
            sym = trim(symbol) || '_' || strip(suffix);
        end;
        else sym = symbol;
        run;

    proc sort data=suffix;
        by patRIL matRIL;
        run;

    proc transpose data=suffix out=cgflip2;
        by patRIL matRIL;
        var col1;
        id sym;
        run;

    data SEM.dsrp_sex_det_sbs_symbol;
        set cgflip2;
        drop _name_;
        run;

/* Create Permanent DSRP side-by-side */
    * Get sex det gene symbols merged with isoform information;
    data suffix;
        set SEM.dsrp_sex_det_stack;
        num = index(cgnumber, '_');
        if num > 0 then do;
            suffix = substr(cgnumber, num+1);
            sym = trim(symbol) || '_' || strip(suffix);
        end;
        else sym = symbol;
        keep patRIL matRIL col1 cgnumber sym;
        run;

    proc sort data=suffix;
        by patRIL matRIL cgnumber;
        run;

    proc sort data=SEM.DSRP_stack;
        by patRIL matRIL cgnumber;
        run;

    data merged;
        merge SEM.DSRP_stack (in=in1) suffix (in=in2);
        by patRIL matRIL cgnumber;
        if in1 and not in2 then sym = cgnumber;
        keep patRIL matRIL col1 sym;
        run;

    proc sort data=merged;
        by patRIL matRIL;
        run;

    proc transpose data=merged out=cgflip3;
        by patRIL matRIL;
        var col1;
        id sym;
        run;

    data SEM.dsrp_sbs_symbol;
        set cgflip3;
        drop _name_;
        run;
        
/* Clean-up */
proc datasets nolist;
    delete cgflip;
    delete cgflip2;
    delete cgflip3;
    delete suffix;
    delete dmel;
    delete flip;
    delete freqs;
    delete merged;
    delete mergeno;
    delete mergeno2;
    delete mergeno3;
    delete mergeok;
    delete mergeok2;
    delete mergeokclean;
    delete mergeokclean2;
    delete oops;
    delete sex;
    run;quit;
