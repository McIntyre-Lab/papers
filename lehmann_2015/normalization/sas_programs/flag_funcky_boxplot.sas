/********************************************************************************
* After making all of the normalization boxplots, there are some samples that
* look strange. To look at this in more detail, I have decided to flag these.
********************************************************************************/

/* By looking at log_uq boxplots I decide these looked a little strange */
data funky;
    length line $ 4;
    length mating_status $ 1;
    length rep $ 1;
    input line $ mating_status $ rep $;
    datalines;
    r158 M 3
    r158 V 3
    r181 V 5
    r38 M 1
    r332 M 1
    r332 M 3
    r356 V 1
    r375 V 1
    r377 M 1
    r377 M 3
    r377 V 1
    r377 V 2
    r391 V 2
    r392 V 2
    r476 M 1
    r476 M 2
    r476 M 3
    r476 V 2
    r476 V 3
    r502 M 2
    r584 M 1
    r584 M 3
    r584 V 2
    r584 V 3
    r737 M 3
    r790 M 1
    r790 M 2
    r790 M 3
    r790 V 1
    r857 M 1
    r857 M 2
    r857 V 1
    r857 V 2
    r897 M 2
    r908 V 2
    w35 M 3
    w69 V 3
    w70 M 2
    w70 M 3
    w74 V 1
    w114 V 2
    ;
    run;

/* Merge onto the complete list of reps and make flag_funky */
proc sort data=funky;
    by line mating_status rep;
    run;

proc sort data=CEGS.combined_design_by_rep;
    by line mating_status rep;
    run;

data CEGS.flag_funky_boxplot;
    merge CEGS.combined_design_by_rep (in=in1) funky (in=in2);
    by line mating_status rep;
    if in2 then flag_funky = 1; else flag_funky=0;
    run;

/* Clean up */
proc datasets nolist;
    delete funky;
    run; quit;
