/*******************************************************************************
* Filename: create_genotype_list.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Create a dataset containing a list of the 68 genotypes for the
* CEGS ASE paper.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/
data CEGS.GENOTYPE_LIST    ;
    infile '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/CEGS_list_68_lines.txt' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=1 ;
    informat line $5. ;
    format line $5. ;
    input
        line $
        ;
    run;

data CEGS.ase_design_file;
    infile '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/CEGS_68_lines_no_tech.txt' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=1 ;
    informat line $5. ;
    informat mating_status $1.;
    informat rep best12. ;
    format line $5. ;
    format mating_status $1.;
    format rep best12. ;
    input
        line $
        mating_status $
        rep
        ;
    run;
