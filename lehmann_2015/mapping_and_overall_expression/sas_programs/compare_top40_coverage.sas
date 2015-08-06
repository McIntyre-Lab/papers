data t40;
    set top40lines;
    flag_t40 = 1;
    drop sum_aln_reads;
    run;

data t40ms;
    set top40lines_ms;
    flag_t40ms = 1;
    drop M V;
    run;

data ht40;
    set hamdi_top40lines;
    flag_ht40 = 1;
    drop sum_apn;
    run;

data ht40ms;
    set hamdi_top40lines_ms;
    flag_ht40ms = 1;
    drop M V;
    run;

proc sort data=t40;
    by line;
    run;
proc sort data=t40ms;
    by line;
    run;
proc sort data=ht40;
    by line;
    run;
proc sort data=ht40ms;
    by line;
    run;

data merged;
    merge t40 t40ms ht40 ht40ms;
    by line;
    run;

data merged2;
    set merged;
    if flag_t40 eq '.' then flag_t40=0;
    if flag_t40ms eq '.' then flag_t40ms=0;
    if flag_ht40 eq '.' then flag_ht40=0;
    if flag_ht40ms eq '.' then flag_ht40ms=0;
    run;

proc export data=merged2 outfile='/home/jfear/mclab/cegs_sergey/reports_external/top40_lines_comparison.csv' dbms=csv replace;
    putnames=yes;
    run;
