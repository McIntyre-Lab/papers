/********************************************************************************
* This script creats a subset of the for analysis table for Lauren and I go sit
* down and go through the modeling.
********************************************************************************/

data FRU.subset_fail_normality;
    retain fusion_id flag_problem_fusion;
    set FRU.foranalysis_no_multi;
     where symbol_cat eq 'Dh31-R1' or symbol_cat eq 'Ir51a' or symbol_cat eq 'Tor' or 
           symbol_cat eq 'dsx' or symbol_cat eq 'fru' or symbol_cat eq 'CG42330';
     if fusion_id eq 'S52794_SI' or fusion_id eq 'S52798_SI' or fusion_id eq 'S47249_SI' or 
        fusion_id eq 'S7129_SI' or fusion_id eq 'S37849_SI' or fusion_id eq 'S38862_SI' or
        fusion_id eq 'S24566_SI'
        then flag_problem_fusion = 1;
     else flag_problem_fusion = 0;
     run;
