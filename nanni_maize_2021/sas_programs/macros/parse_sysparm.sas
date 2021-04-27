/*******************************************************************************
* Filename: parse_sysparam.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Parse multiple sysparm options in the format of:
*   "name1=value1,name2=value2,name3=value3" use with an include not a macro
*   call
*
*******************************************************************************/

data _null_;
    length sysparm express param value $300;
    sysparm = symget('sysparm');
    do i=1 to &cnt;
        express = left(scan(sysparm, i, ','));
        param = left(scan(express, 1, '='));
        value = left(scan(express, 2, '='));
        call symput(param, trim(left(value)));
    end;
    run;
