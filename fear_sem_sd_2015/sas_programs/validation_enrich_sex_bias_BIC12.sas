/*******************************************************************************
* Filename: validation_enrich_sex_bias_BIC12.sas
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
    set SEM.validation_set_BIC12;
    drop model modelnum path BIC;
    run;

proc sort data=valid nodupkey out=genelevel;
    by primary_fbgn;
    run;

proc freq data=genelevel;
    tables flag_grn_expansion;
    run; * 754 / 8810 = 8.56  %;

proc freq data=genelevel;
    tables flag_grn_expansion*flag_mcintyre_sex_bias_splice;
    run; * 11 / 38 = 28.95  %;

/* Run two sample t-test to compare proporitions */
/*

Ho: P1 => P2
Ha: P1 > P2
P1 = 11/38                                     = 0.2895
P2 = 754/8810                                  = 0.0856
Pt = (11 + 754) / (38 + 8810)                  = 0.0864
SE = sqrt(Pt * (1 - Pt) * (1/8810 + 1/38))     = 0.0457
z = (P1 - P2) / SE                             = 4.4625
P(Z > z) = pnorm(z, lower.tail=FALSE)          = <0.0001

*/
