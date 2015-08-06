/********************************************************************************
* There was some discussion to create a list of the top 40 lines based on gene
* expression. So this set of programs identifies lines with the most coverage
********************************************************************************/

libname cegs '!MCLAB/cegs_sergey/sas_data';
libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
libname datadump '!SASLOC2/cegs_sergey/sasdata';
libname mixed '!SASLOC2/cegs_sergey/mixed_sasdata';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Identify lines with the most coverage based on all genes */
    * This script imports the alignment summary and calculates the top 40 lines
    * based on total coverage and based on MV coverage separately.
    *
    * INFILE: !MCLAB/cegs_sergey/pipeline_output/aln_summary_fb551_non-redundant_fusions_20130912.csv
    *
    * DATASET: WORK.top40lines
    *          WORK.top40lines_ms
    *
    * FILES: !MCLAB/cegs_sergey/reports_external/top_40_lines_coverage.csv
    *        !MCLAB/cegs_sergey/reports_external/top_40_lines_w_mating_status_coverage.csv
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/best_coverage_lines.sas';

/* Compare different top 40 lists */
    * Since I have several different top 40 lists, I wanted to compare things
    * and see how similar they were.
    *
    * INPUT: WORK.top40lines
    *        WORK.top40lines_ms
    *        WORK.hamdi_top40lines
    *        WORK.hamdi_top40lines_ms
    * FILES: !MCLAB/cegs_sergey/reports_external/top40_lines_comparison.csv
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/compare_top40_coverage.sas';
