/********************************************************************************
* This script imports tables that Michelle supplied. These tables are genes
* lists of genes that are repressed or induced by Fru. There are 6 tables
* representing Fru A or B or C
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';

%macro import_files(input_file,out_file,guessing);
proc import datafile = &input_file
            out=&out_file
            dbms=TAB replace;
            getnames=no;
            guessingrows=&guessing;
            run; *602 obs same as tsv file;
%mend import_files;

%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/FruA male induced.tab',FruA_ind,1897);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru A male repressed.tab',FruA_rep,268);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru B male induced.tab',FruB_ind,1653);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru B male repressed.tab',FruB_rep,343);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/FruC male induced.tab',FruC_ind,2158);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru C male repressed.tab',FruC_rep,356);
