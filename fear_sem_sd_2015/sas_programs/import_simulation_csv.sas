
/* Split sysparm */
data _null_;
    length sysparm express param value $200;
    sysparm = symget('sysparm');
    do i=1 to 1;
        express = left(scan(sysparm, i, ','));
        param = left(scan(express, 1, '='));
        value = left(scan(express, 2, '='));
        call symput(param, trim(left(value)));
    end;
    run;

/* Combine models */
libname TMP "&lib1";

data TMP.simulation;
    set rdata;
    run;
