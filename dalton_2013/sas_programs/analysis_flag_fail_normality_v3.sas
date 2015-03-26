libname fru '!MCLAB/Fru_network/sasdata';

ods listing close;
proc univariate data=fru.resid_by_symbol normal plots;
    by symbol_cat;
    var resid;
    ods output TestsForNormality=normtest;
    run;
ods listing;

data fru.flag_fail_normality_by_symbol;
    set normtest;
    where test="Shapiro-Wilk";
    if pValue eq . then flag_fail_normality=.;
    else if pValue <= .05 then flag_fail_normality=1;
    else flag_fail_normality=0;
    keep symbol_cat flag_fail_normality;
    run;

proc freq data=fru.flag_fail_normality_by_symbol;
    tables flag_fail_normality;
    run;
    

/* resid checks */

/* normal */
proc gplot data = fru.resid_by_symbol (where=(symbol_cat='14-3-3epsilon')) ;
    plot resid*pred ;
run;
quit;

proc gplot data = fru.resid_by_symbol (where=(symbol_cat='26-29-p')) ;
    plot resid*pred ;
run;
quit;

/* fail normallity */
proc gplot data = fru.resid_by_symbol (where=(symbol_cat='128up')) ;
    plot resid*pred ;
run;
quit;

proc gplot data = fru.resid_by_symbol (where=(symbol_cat='msn')) ;
    plot resid*pred ;
run;
quit;
