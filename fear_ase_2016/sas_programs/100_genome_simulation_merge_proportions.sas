/*******************************************************************************
* Filename: 100_genome_simulation_merge_proportions.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Take the proporation of tester / total from each genotype and
* merge into a sbs.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname MYCEGS '!HOME/storage/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

proc sort data=all_ase;
    by fusion_id;
    run;

data wide;
    set all_ase;
    by fusion_id;
    keep fusion_id line1-line100;
    retain line1-line100;
    array aprop(1:100) line1-line100;
    if first.fusion_id then do;
        do i = 1 to 100;
            aprop(i) = .;
        end;
    end;
    num = substr(line, 5);
    aprop( num ) = prop_tester;
    if last.fusion_id then output;
    run;

proc export data=wide outfile='!MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/fb551_100_genome_bias_wide.csv' dbms=csv replace;
    putnames=yes;
    run;

data CEGS.fb551_100_genome_flag_line_bias;
    set all_ase;
    by fusion_id;
    keep fusion_id line1-line100;
    retain line1-line100;
    array aprop(1:100) line1-line100;
    if first.fusion_id then do;
        do i = 1 to 100;
            aprop(i) = 0;
        end;
    end;
    num = substr(line, 5);
    if prop_tester < 0.5 then flag = -1;  * towards line;
    else if prop_tester > 0.5 then flag = 1; * towards tester;
    else flag = 0;
    aprop( num ) = flag;
    if last.fusion_id then output;
    run;

proc export data=CEGS.fb551_100_genome_flag_line_bias outfile='!MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/fb551_100_genome_flag_line_bias.csv' dbms=csv replace;
    putnames=yes;
    run;

proc datasets nolist;
    delete all_ase;
    run; quit;
