/*
    libname CEGS '!MCLAB/cegs_sergey/sas_data';
    libname CEGLOCAL '!SASLOC1/cegs_sergey/sasdata';
    libname DMEL '!MCLAB/useful_dmel_data/flybase551/sasdata';
    filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
    options SASAUTOS=(sasautos mymacros);
*/

/* plate k */
data platek;
    set CEGS.platek_ccfus_stack;
    rename mapped_reads = platek_mapped_reads;
    keep line mating_status rep mapped_reads;
    run;

proc sort data=platek nodupkey;
    by _all_;
    run;

* Sum mapped reads;
proc means data=platek noprint;
    by line mating_status;
    output out=sums sum(platek_mapped_reads)=platek_total_mapped;
    run;

data sums2;
    set sums;
    drop _type_ _freq_;
    run;

* count how many reps;
proc freq data=platek noprint;
    table line*mating_status /out=freq;
    run;

data freq2;
    set freq;
    rename count=platek_num_rep;
    label count = ' ';
    drop percent;
    run;

* Put it all together;
proc sort data=freq2;
    by line mating_status;
    run;

proc sort data=platek;
    by line mating_status;
    run;

proc sort data=sums2;
    by line mating_status;
    run;

data platek_w_cnt;
    merge platek (in=in1) freq2 (in=in2) sums2 (in=in3);
    by line mating_status;
    run;


/* before plate k */
data prior;
    set CEGS.ccfus_stack;
    keep line mating_status rep mapped_reads;
    run;

proc sort data=prior nodupkey;
    by _all_;
    run;

* Sum mapped reads;
proc means data=prior noprint;
    by line mating_status;
    output out=sums sum(mapped_reads)=total_mapped;
    run;

data sums2;
    set sums;
    drop _type_ _freq_;
    run;

* count how many reps;
proc freq data=prior noprint;
    table line*mating_status /out=freq;
    run;

data freq2;
    set freq;
    rename count=num_rep;
    label count = ' ';
    drop percent;
    run;

* Put it all together;
proc sort data=freq2;
    by line mating_status;
    run;

proc sort data=prior;
    by line mating_status;
    run;

proc sort data=sums2;
    by line mating_status;
    run;

data prior_w_cnt;
    merge prior (in=in1) freq2 (in=in2) sums2(in=in3);
    by line mating_status;
    run;


/* Merge everything together */
* rep level;
proc sort data=platek_w_cnt;
    by line mating_status rep;
    run;

proc sort data=prior_w_cnt;
    by line mating_status rep;
    run;

data merged;
    merge platek_w_cnt (in=in1) prior_w_cnt (in=in2);
    by line mating_status rep;
    if in2 and not in1 then do;
        platek_mapped_reads = 0;
    end;
    drop platek_total_mapped platek_num_rep total_mapped num_rep;
    run;

data CEGS.platek_comparison_rep;
    set merged;
    run;

* overall level;
data pk;
    set platek_w_cnt;
    drop rep platek_mapped_reads;
    run;

proc sort data=pk nodupkey;
    by line mating_status;
    run;

data pr;
    set prior_w_cnt;
    drop rep mapped_reads;
    run;

proc sort data=pr nodupkey;
    by line mating_status;
    run;

data merged_tot;
    merge pk (in=in1) pr (in=in2);
    by line mating_status;
    if in2 and not in1 then do;
        platek_num_rep = 0;
        platek_total_mapped = 0;
    end;
    if platek_total_mapped gt total_mapped then flag_platek_gt = 1;
    else flag_platek_gt = 0;
    run;

data CEGS.platek_comparison_tot;
    set merged_tot;
    run;

/* Export files */
proc export data=CEGS.platek_comparison_rep outfile='!MCLAB/cegs_sergey/reports/platek_comparison_rep.csv' dbms=csv replace;
putname=yes;
run;

proc export data=CEGS.platek_comparison_tot outfile='!MCLAB/cegs_sergey/reports/platek_comparison_total.csv' dbms=csv replace;
putname=yes;
run;


/* Clean up */
proc datasets nolist;
    delete freq;
    delete freq2;
    delete merged;
    delete merged_tot;
    delete pk;
    delete platek;
    delete platek_w_cnt;
    delete pr;
    delete prior;
    delete prior_w_cnt;
    delete sums;
    delete sums2;
    run; quit;


