/*******************************************************************************
* Filename: cis_eq_for_clustering.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create Split up data matrices and Export */
    proc sort data=CEGS.cis_eq;
        by fusion_id;
        run;

    * Cis-line Mated;
    proc transpose data=CEGS.cis_eq out=cism;
        where mating_status eq 'M';
        by fusion_id;
        var cis_line;
        id line;
        run;

    data CEGS.cis_line_effect_mated;
        set cism;
        drop _name_;
        run;

    proc export data=CEGS.cis_line_effect_mated outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/cis_line_effect_mated.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Cis-line Virgin;
    proc transpose data=CEGS.cis_eq out=cisv;
        where mating_status eq 'V';
        by fusion_id;
        var cis_line;
        id line;
        run;

    data CEGS.cis_line_effect_virgin;
        set cisv;
        drop _name_;
        run;

    proc export data=CEGS.cis_line_effect_virgin outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/cis_line_effect_virgin.csv' dbms=csv replace;
        putnames=yes;
        run;

    * trans-line Mated;
    proc transpose data=CEGS.cis_eq out=transm;
        where mating_status eq 'M';
        by fusion_id;
        var trans_line;
        id line;
        run;

    data CEGS.trans_line_effect_mated;
        set transm;
        drop _name_;
        run;

    proc export data=CEGS.trans_line_effect_mated outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/trans_line_effect_mated.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Cis-line Virgin;
    proc transpose data=CEGS.cis_eq out=transv;
        where mating_status eq 'V';
        by fusion_id;
        var trans_line;
        id line;
        run;

    data CEGS.trans_line_effect_virgin;
        set transv;
        drop _name_;
        run;

    proc export data=CEGS.trans_line_effect_virgin outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/trans_line_effect_virgin.csv' dbms=csv replace;
        putnames=yes;
        run;


/* Removing Missing and output for MMC */
    * Delete rows with missing values;
    data mclean;
        set CEGS.cis_line_effect_mated;
        if nmiss(of _numeric_) >0 then delete;
        run;

    * Transpose matrix so that genotypes are rows;
    proc transpose data=mclean out=mmmc;
        var :;
        id fusion_id;
        run;

    proc export data=mmmc outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/mmc/cis_line_effect_for_mmc_mated.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Delete rows with missing values;
    data vclean;
        set CEGS.cis_line_effect_virgin;
        if nmiss(of _numeric_) >0 then delete;
        run;

    * Transpose matrix so that genotypes are rows;
    proc transpose data=vclean out=vmmc;
        var :;
        id fusion_id;
        run;

    proc export data=vmmc outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/mmc/cis_line_effect_for_mmc_virgin.csv' dbms=csv replace;
        putnames=yes;
        run;
