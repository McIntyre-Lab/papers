/*******************************************************************************
* Filename: 100_genome_simulation_create_drop_list.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: I made a drop list in python, I now want to import the flags into
* SAS. See $PROJ/scripts/100_genome_simulation/genome_ambiguity_summary.ipynb
* for more information.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/


proc import datafile='!MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/flag_exonic_region_w_and_wo_bias_100_genome_simulation.xls' out=CEGS.flag_exonic_regions_100_genome dbms=xls replace;
    getnames=yes;
    run;

data CEGS.exon_drop_list_100_genome;
    set CEGS.flag_exonic_regions_100_genome;
    where flag_exons_bias_in_all_simulated eq 1;
    keep fusion_id;
    run;
