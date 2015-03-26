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

%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru A female induced.tab',FruA_ind,162);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru A female repressed.tab',FruA_rep,265);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru B female induced.tab',FruB_ind,167);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru B female repressed.tab',FruB_rep,302);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru C female induced.tab',FruC_ind,289);
%import_files('!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Fru C female repressed.tab',FruC_rep,244);
