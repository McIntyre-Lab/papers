/* Split sysparm */
data _null_;
    length sysparm express param value $200;
    sysparm = symget('sysparm');
    do i=1 to 3;
        express = left(scan(sysparm, i, ','));
        param = left(scan(express, 1, '='));
        value = left(scan(express, 2, '='));
        call symput(param, trim(left(value)));
    end;
    run;

/* Combine models */
libname TMP "&lib1";
libname STORE "&lib2";

data STORE.&gene;
    set TMP.gene_:;
    label BIC = ' ';
    run;

proc sort data=STORE.&gene;
    by BIC;
    run;
