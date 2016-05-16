/*******************************************************************************
* Filename: ase_summarize_eqtl_tf.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Sergey has given me a list of genes that are affected by the
* eQTLs with snps in TF binding sites. I want to merge to this list and see if
* any of these have AI.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Import eQTL TF list */
    proc import datafile='!MCLAB/cegs_ase_paper/external_data/sergey_eqtl_transcription_factor_list.txt' out=eqtl dbms=csv replace;
        getnames=yes;
        run;

/* Select exonic regions with AI */ 
    data ai;
        set CEGS.clean_ase_sbs;
        if flag_AI_combined_m eq 1 or flag_AI_combined_v eq 1;
        run;

/* Merge FBGN to exonic region ids */
    proc sort data=DMEL551.fb551_si_fusion_2_gene_id;
        by fusion_id;
        run;

    proc sort data=ai;
        by fusion_id;
        run;

    data aimerged;
        merge ai (in=in1) DMEL551.fb551_si_fusion_2_gene_id (in=in2 keep=fusion_id FBgn);
        by fusion_id;
        if in1;
        run;

/* Merge eQTL to ai */
    proc sort data=eqtl;
        by primary_fbgn;
        run;

    proc sort data=aimerged;
        by fbgn;
        run;

    data qtlmerged;
        merge aimerged (in=in1 rename=(fbgn=primary_fbgn)) eqtl (in=in2 );
        by primary_fbgn;
        if in1 and in2;
        run;

/* Look at counts */
proc sort data=qtlmerged;
    by primary_fbgn;
    run;

proc freq data=qtlmerged ;
    by primary_fbgn;
    tables flag_AI_combined_m;
    tables flag_AI_combined_v;
    run;
