/*******************************************************************************
* Filename: ase_summarize_ase_discordance_mating_status.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: We are interested in exonic regions/genes that switch AI
* direction with the environment. This could be something of biological
* interest. Pulls out these exonic regions and does some basic counts.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
*/

data ai; * exonic regions with AI;
    set CEGS.clean_ase_sbs;
    if (flag_AI_combined_m eq 1 and flag_AI_combined_v eq 1);
    drop flag_AI_combined_m flag_AI_combined_v;
    run;

data dis; * discordant exonic regions;
    set ai;
    if (q5_mean_theta_m < 0.5 and q5_mean_theta_v > 0.5) or
       (q5_mean_theta_m > 0.5 and q5_mean_theta_v < 0.5) then flag_AI_and_discord = 1;
    else flag_AI_and_discord = 0;
    distance = abs(q5_mean_theta_m - q5_mean_theta_v);
    run;

/* Table of fusions with discord */
    proc sort data=dis;
        by fusion_id;
        run;

    proc freq data=dis noprint;
        by fusion_id;
        tables flag_AI_and_discord /out=genoDis; * counts of level of discordance;
        run;

    proc transpose data=genoDis out=genoDisTrans prefix=AI_discord_ suffix=_cnt;
        by fusion_id;
        var count;
        id flag_AI_and_discord;
        run;

    data genoDisTrans;
        set genoDisTrans;
        if AI_discord_0_cnt eq . then AI_discord_0_cnt = 0;
        if AI_discord_1_cnt eq . then AI_discord_1_cnt = 0;
        keep fusion_id AI_discord_0_cnt AI_discord_1_cnt;
        run;

    * Are there any fusions with many genotypes showing discord;
    proc freq data=genoDisTrans;
        tables AI_discord_1_cnt;
        run;


/* Table of genotypes with discord */
    proc sort data=dis;
        by line;
        run;

    proc freq data=dis noprint;
        by line;
        tables flag_AI_and_discord /out=fusDis; * counts of level of discordance;
        run;

    proc transpose data=fusDis out=fusDisTrans prefix=AI_discord_ suffix=_cnt;
        by line;
        var count;
        id flag_AI_and_discord;
        run;

    data fusDisTrans;
        set fusDisTrans;
        if AI_discord_0_cnt eq . then AI_discord_0_cnt = 0;
        if AI_discord_1_cnt eq . then AI_discord_1_cnt = 0;
        keep line AI_discord_0_cnt AI_discord_1_cnt;
        run;

    * Are there any genotypes with many fusions showing discord;
    proc freq data=fusDisTrans;
        tables AI_discord_1_cnt;
        run;


/* Export Table */
    proc export data=dis outfile='!MCLAB/cegs_ase_paper/pipeline_output/ase_summary/environment_discordance.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
    proc datasets nolist;
        delete AI;
        delete dis;
        delete genodis;
        delete genodisTrans;
        delete fusdis;
        delete fusdisTrans;
        run; quit;

