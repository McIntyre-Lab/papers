/********************************************************************************
* Collapse all isoforms to the gene level, so that I can compare it to
* the cegs dataset.
********************************************************************************/

/* Transpose data */
proc transpose data=SEM.dsrp_sbs_symbol out=flip;
    by patRIL matRIL;
    var : ;
    run;


/* Remove isoform info from gene symbol */
    * Remember that fl_2_d has multiple '_' so treat it separately;

    data sym;
        retain patRIL matRIL symbol _name_ col1;
        length symbol $11.;
        set flip;
        if prxmatch('/fl_2_d/', _name_) then symbol = 'fl_2_d';
        else do;
            num = index(_name_, '_');
            if num > 0 then symbol = substr(_name_ , 1, num-1);
            else symbol = _name_;
        end;
        drop num;
        run;

/* Average Across patRIL matRIL symbol */
    proc sort data=sym;
        by patRIL matRIL symbol;
        run;

    proc means data=sym noprint;
        by patRIL matRIL symbol;
        output out=means mean(col1)=means;
        run;

    data SEM.dsrp_stack_gene_level_sym;
        set means;
        rename means=exp;
        drop _type_ _freq_;
        run;

/* Flip back to SBS */
    proc transpose data=means out=flipback;
        by patRIL matRIL;
        var means;
        id symbol;
        run;

/* Make permanent dataset */
    data SEM.dsrp_sbs_gene_level_sym;
        set flipback;
        drop _name_;
        run;

    /* Use JMP to export dataset, too many columns
    proc export data=SEM.dsrp_sbs_gene_level_sym outfile='!MCLAB/cegs_sem_sd_paper/exported_data/dspr_by_gene_sbs.csv' dbms=csv replace;
        putnames=yes;
        run;
    */

/* Clean up */
    proc datasets nolist;
        delete flip;
        delete flipback;
        delete means;
        delete sym;
        run; quit;
