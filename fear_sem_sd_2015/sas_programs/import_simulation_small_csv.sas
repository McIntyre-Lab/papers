/* Split sysparm */
data _null_;
    length sysparm express param value $200;
    sysparm = symget('sysparm');
    do i=1 to 2;
        express = left(scan(sysparm, i, ','));
        param = left(scan(express, 1, '='));
        value = left(scan(express, 2, '='));
        call symput(param, trim(left(value)));
    end;
    run;

/* Combine models */
libname TMP "&lib";

proc import datafile="&mydat" out=TMP.simulated_data dbms=csv replace;
    getnames=yes;
    run;
