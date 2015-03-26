/*
 * This script is used to merege the needed information together to make a
 * reference for use in the r_wiggle script.
 */

%let MCLAB=!MCLAB

libname dmel530 '&MCLAB/useful_dmel_data/Flybase 5.30/dmel_annotation';

proc import out=work.fusions
            datafile='&MCLAB/useful_dmel_data/Flybase 5.30/dmel_annotation'
            dbms=csv replace;
            datarows=2;
            getnames=yes;
            run;
