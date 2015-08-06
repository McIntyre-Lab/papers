/*******************************************************************************
* Filename: validation_enrich_sex_bias.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Sex det pathway is a splicing cascade, so we may expect that
* genes previously shown as sex bias, may be effected.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/
 
data valid;
    set SEM.validation_set;
    drop model modelnum path BIC;
    run;

proc sort data=valid nodupkey out=genelevel;
    by primary_fbgn;
    run;

proc freq data=genelevel;
    tables flag_grn_expansion;
    run; * 1390 / 8810 = 15.78 %;

proc freq data=genelevel;
    tables flag_grn_expansion*flag_mcintyre_sex_bias_splice;
    run; * 18 / 38 = 47.36 %;

/* Run two sample t-test to compare proporitions */
/*
Ho: P1 = P2
Ha: P1 ne P2
P1 = 18/38 = .4736
P2 = 1390/8810 = .1578
Pt = (18 + 1390) / (38 + 8810) = .1591

z = (.4736 - .1578 - 0)/ sqrt(.1591(1-.1591)(1/38 + 1/8810)) = .3158 / sqrt(.1338(.0264)) = 5.31

p-value = <0.0001

*/
