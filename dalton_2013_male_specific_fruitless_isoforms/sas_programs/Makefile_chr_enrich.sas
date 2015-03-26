/********************************************************************************
* This is the Makefile for fru chromosomal Analysis
********************************************************************************/

libname fru '!MCLAB/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* DATA PREP */
    /** MALE **/
        /*** Import Michelles List ***/
            * Import Gene list provided by Michelle. These genes are considered to be
            * significantly differentially expressed.
            *
            * INPUT:
            *   /mclab/Fru_network/original_data/FruA%20male%20induced.tab
            *   /mclab/Fru_network/original_data/Fru%20A%20male%20repressed.tab
            *   /mclab/Fru_network/original_data/Fru%20B%20male%20induced.tab
            *   /mclab/Fru_network/original_data/Fru%20B%20male%20repressed.tab
            *   /mclab/Fru_network/original_data/FruC%20male%20induced.tab
            *   /mclab/Fru_network/original_data/Fru%20C%20male%20repressed.tab"
            *
            * DATASETS:
            *   WORK.FruA_ind
            *   WORK.FruA_rep
            *   WORK.FruB_ind
            *   WORK.FruB_rep
            *   WORK.FruC_ind
            *   WORK.FruC_rep
            *
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_import_repressed_induced_files_male.sas';

        /*** Create Union ***/
            * Create a union list of repressed genes and a union list of induced genes.
            *
            * INPUT:
            *   WORK.FruA_ind
            *   WORK.FruA_rep
            *   WORK.FruB_ind
            *   WORK.FruB_rep
            *   WORK.FruC_ind
            *   WORK.FruC_rep
            *
            * DATASETS:
            *   WORK.FruA_ind_nodups = 752 obs
            *   WORK.FruA_rep_nodups = 204 obs
            *   WORK.FruB_ind_nodups = 739 obs
            *   WORK.FruB_rep_nodups = 259 obs
            *   WORK.FruC_ind_nodups = 927 obs
            *   WORK.FruC_rep_nodups = 295 obs
            *
            *   FRU.male_induced_union   = 1217 obs
            *   FRU.male_repressed_union = 554 obs
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_create_union_male_repressed_induced.sas';

        /*** Flag for Induced and Repressed ***/
            * Create a union list of repressed genes and a union list of induced genes.
            *
            * INPUT:
            *   FRU.male_induced_union   = 1217 obs
            *   FRU.male_repressed_union = 554 obs
            *
            * DATASETS:
            *   WORK.all_merged
            *   WORK.oops
            *   FRU.flag_x_induced_repressed_male
            *   FRU.flag_ind_rep_w_het_male
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_flag_repressed_induced_male_v4.sas';

        /*** Export flag_x_induced_repressed_male ***/
            * Export the flags to a CSV
            *
            * INPUT: FRU.flag_x_induced_repressed_male
            *
            * OUTPUT: '!MCLAB/arbeitman_fru_network/reports_internal/flag_x_induced_repressed_male.csv'
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_flag_x_induced_repressed_male.sas';

    /** FEMALE **/
        /*** Import Michelles List ***/
            * Import Gene list provided by Michelle. These genes are considered to be
            * significantly differentially expressed.
            *
            * INPUT:
            *   /mclab/Fru_network/original_data/FruA%20female%20induced.tab
            *   /mclab/Fru_network/original_data/Fru%20A%20female%20repressed.tab
            *   /mclab/Fru_network/original_data/Fru%20B%20female%20induced.tab
            *   /mclab/Fru_network/original_data/Fru%20B%20female%20repressed.tab
            *   /mclab/Fru_network/original_data/FruC%20female%20induced.tab
            *   /mclab/Fru_network/original_data/Fru%20C%20female%20repressed.tab"
            *
            * DATASETS:
            *   WORK.FruA_ind
            *   WORK.FruA_rep
            *   WORK.FruB_ind
            *   WORK.FruB_rep
            *   WORK.FruC_ind
            *   WORK.FruC_rep
            *
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_import_repressed_induced_files_female.sas';

        /*** Create Union ***/
            * Create a union list of repressed genes and a union list of induced genes.
            *
            * INPUT:
            *   WORK.FruA_ind
            *   WORK.FruA_rep
            *   WORK.FruB_ind
            *   WORK.FruB_rep
            *   WORK.FruC_ind
            *   WORK.FruC_rep
            *
            * DATASETS:
            *   WORK.FruA_ind_nodups = 111 obs
            *   WORK.FruA_rep_nodups = 183 obs
            *   WORK.FruB_ind_nodups = 117 obs
            *   WORK.FruB_rep_nodups = 237 obs
            *   WORK.FruC_ind_nodups = 167 obs
            *   WORK.FruC_rep_nodups = 198 obs
            *
            *   FRU.female_induced_union   = 267 obs
            *   FRU.female_repressed_union = 462 obs
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_create_union_female_repressed_induced.sas';

        /*** Flag for Induced and Repressed ***/
            * Create a union list of repressed genes and a union list of induced genes.
            *
            * INPUT:
            *   FRU.female_induced_union   = 1217 obs
            *   FRU.female_repressed_union = 554 obs
            *
            * DATASETS:
            *   WORK.all_merged
            *   WORK.oops
            *   FRU.flag_x_induced_repressed_female
            *   FRU.flag_ind_rep_w_het_female
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_flag_repressed_induced_female_v4.sas';

        /*** Export flag_x_induced_repressed_female ***/
            * Export the flags to a CSV
            *
            * INPUT: FRU.flag_x_induced_repressed_female
            *
            * OUTPUT: '!MCLAB/arbeitman_fru_network/reports_internal/flag_x_induced_repressed_female.csv'
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_flag_x_induced_repressed_female.sas';

    /** NULL MALE **/
        /*** Import MY List ***/
            * Import Gene list I created above. These genes are considered to be
            * significantly differentially expressed.
            *
            * INPUT:
            *   /mclab/arbeitman_fru_network/exported_data_from_michelle/Induced_Fru_m_null_jmf.tab
            *   /mclab/arbeitman_fru_network/exported_data_from_michelle/Repressed_Fru_m_null_jmf.tab
            *
            * DATASETS:
            *   WORK.null_ind_nodups
            *   WORK.null_rep_nodups
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_import_repressed_induced_files_null_male_v2.sas';

        /*** Create Union and Flag Induced and Repressed ***/
            * Create a union list of repressed genes and a union list of induced genes.
            *
            * INPUT:
            *   WORK.null_ind_nodups
            *   WORK.null_rep_nodups
            *
            * DATASETS: FRU.flag_x_ind_rep_null_male
            *
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_flag_repressed_induced_null_male_v2.sas';

        /*** Export flag_x_induced_repressed_null_male ***/
            * Export the flags to a CSV
            *
            * INPUT: FRU.flag_x_ind_rep_null_male
            *
            * OUTPUT: '!MCLAB/arbeitman_fru_network/reports_internal/flag_x_induced_repressed_female.csv'
            ;
            %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_flag_x_induced_repressed_null_male.sas';

/* Create Master Flags for Induced Repressed */
    * I want to create a master flag list of
    * {male,female,null},{fruA,fruB,fruC},{induced,repressed,nodiff}
    *
    * INPUT: FRU.Flag_x_induced_repressed_male;
    *        FRU.Flag_x_induced_repressed_female;
    *        FRU.Flag_x_ind_rep_null_male;
    *
    * DATASET: FRU.flag_ind_rep
    ;
    %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_create_master_flags_ind_rep.sas';









/* Calculate Chromosomal Enrichment Overexpression MALE and FEMALE */
    /** Combine Euchromatin and Heterocromatin **/
        * This script performs a fisher's exact test for x chromosome enrichment;
        *
        * INPUT:
        *   FRU.flag_x_induced_repressed_male
        *   FRU.flag_x_induced_repressed_female
        *
        * OUTPUT:
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tables.csv
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tests.csv
        ;
        %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_tests_v2.sas';

    /** Separate Euchromatin and Heterocromatin **/
        * This script performs a fisher's exact test for x chromosome enrichment;
        *
        * INPUT:
        *   FRU.flag_ind_rep_w_het_male
        *   FRU.flag_ind_rep_w_het_female
        *
        * OUTPUT:
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tables_w_het.csv
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tests_w_het.csv
        ;
        %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_tests_w_het.sas';

    /** Split up different Fru's Combine Euchromatin and Heterocromatin **/
        * This script performs a fisher's exact test for x chromosome enrichment;
        *
        * INPUT:
        *   FRU.flag_x_induced_repressed_male
        *   FRU.flag_x_induced_repressed_female
        *
        * OUTPUT:
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tables.csv
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tests.csv
        ;
        %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_tests_by_fru.sas';

    /** Split up different Fru's Separate Euchromatin and Heterocromatin **/
        * This script performs a fisher's exact test for x chromosome enrichment;
        *
        * INPUT:
        *   FRU.flag_ind_rep_w_het_male
        *   FRU.flag_ind_rep_w_het_female
        *
        * DATASET:
        *   FRU.flag_x_ind_rep_null_male
        *
        * OUTPUT:
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tables_w_het.csv
        *   /mclab/arbeitman_fru_network/reports/chrom_enrichment_tests_w_het.csv
        ;
        %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_tests_w_het_by_fru.sas';

/* Calculate Chromosomal Enrichment Null MALE */
    * 
    * INPUT: 
    *   FRU.flag_x_ind_rep_null_male
    *
    * OUTPUT:
    *   ~/mclab/arbeitman_fru_network/reports_external/chrom_enrichment_tables_null_v2.csv
    *   ~/mclab/arbeitman_fru_network/reports_external/chrom_enrichment_tests_null_v2.csv
    ;
    %include '!MCLAB/arbeitman_fru_network/sas_programs/chr_enrich_tests_null.sas';

/* Compare Overexpression vs Null */
    * 
    ;

data ov_ind;
    set FRU.flag_x_induced_repressed_male ;
    rename flag_induced = ov_ind;
    rename flag_repressed = ov_rep;
    keep primary_fbgn symbol flag_induced flag_repressed;
    run;

data null_ind;
    set FRU.flag_x_ind_rep_null_male ;
    rename flag_induced = null_ind;
    rename flag_repressed = null_rep;
    keep primary_fbgn symbol flag_induced flag_repressed;
    run;

proc sort data=ov_ind;
    by primary_fbgn;
    run;

proc sort data=null_ind;
    by primary_fbgn;
    run;

data over_vs_null;
    merge ov_ind (in=in1) null_ind (in=in2);
    by primary_fbgn;
    run;

proc freq data=over_vs_null;
    tables ov_ind*null_ind /agree chisq expected ;
    exact fisher chisq;
    run;

