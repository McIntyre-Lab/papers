/* Calculate Ranks on the Entire DSRP dataset */
    proc sort data=SEM.dsrp_stack;
        by patRIL matRIL;
        run;
        
    proc rank data=SEM.dsrp_stack out=SEM.dsrp_stack_rank;
        by patRIL matRIL;
        var col1;
        ranks expRank;
        run;

