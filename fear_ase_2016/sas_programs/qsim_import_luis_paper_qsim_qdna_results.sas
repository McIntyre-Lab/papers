/*******************************************************************************
* Filename: qsim_import_luis_paper_qsim_qdna_results.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import the results from the Luis paper.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/


data WORK.LUIS    ;
    infile '!MCLAB/cegs_ase_paper/pipeline_output/qsim/luis_paper/luis_nofilter_SNPout_NBandPG_Modelsonlyp_simnotahalfwithmeanssmltd_freq_mel_allfusions.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat gene_symbol_cat $129. ;
        informat FUSION_ID $9. ;
        informat RNA_APN_detected best32. ;
        informat S_RNA best32. ;
        informat TotalRNA best32. ;
        informat PER_S_RNA best32. ;
        informat S_DNA best32. ;
        informat TotalDNA best32. ;
        informat PER_S_DNA best32. ;
        informat p_sim best32. ;
        informat mean_CI best32. ;
        informat q025_CI best32. ;
        informat q975_CI best32. ;
        informat mean_PG_q_DNA best32. ;
        informat q025_PG_q_DNA best32. ;
        informat q975_PG_q_DNA best32. ;
        informat mean_PG_q_sim best32. ;
        informat q025_PG_q_sim best32. ;
        informat q975_PG_q_sim best32. ;
        informat mean__PG_q_ahalf best32. ;
        informat q025_PG_q_ahalf best32 ;
        informat q975_PG_q_ahalf best32. ;
        format gene_symbol_cat $129. ;
        format FUSION_ID $9. ;
        format RNA_APN_detected best12. ;
        format S_RNA best12. ;
        format TotalRNA best12. ;
        format PER_S_RNA best12. ;
        format S_DNA best12. ;
        format TotalDNA best12. ;
        format PER_S_DNA best12 ;
        format p_sim best12. ;
        format mean_CI best12. ;
        format q025_CI best12. ;
        format q975_CI best12. ;
        format mean_PG_q_DNA best12. ;
        format q025_PG_q_DNA best12. ;
        format q975_PG_q_DNA best12. ;
        format mean_PG_q_sim best12. ;
        format q025_PG_q_sim best12. ;
        format q975_PG_q_sim best12. ;
        format mean__PG_q_ahalf best12. ;
        format q025_PG_q_ahalf best12. ;
        format q975_PG_q_ahalf best12. ;
        input
            gene_symbol_cat $
            FUSION_ID $
            RNA_APN_detected
            S_RNA
            TotalRNA
            PER_S_RNA 
            S_DNA
            TotalDNA
            PER_S_DNA 
            p_sim 
            mean_CI 
            q025_CI 
            q975_CI 
            mean_PG_q_DNA 
            q025_PG_q_DNA 
            q975_PG_q_DNA 
            mean_PG_q_sim 
            q025_PG_q_sim 
            q975_PG_q_sim 
            mean__PG_q_ahalf 
            q025_PG_q_ahalf 
            q975_PG_q_ahalf
            ;
        run;

data CEGS.luis_pg_results;
    set luis;
    run;
