/********************************************************************************
* This script imports tables that Michelle supplied. These tables are genes
* lists of genes that are repressed or induced by Fru. There are 6 tables
* representing Fru A or B or C
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';

%macro import_files(input_file,out_file,guessing);

    data &out_file    ;
        infile &input_file delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat FBgn $11. ;
        format FBgn $11. ;
        input FBgn $ ;
        run;

%mend import_files;

%import_files(input_file='!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Induced_Fru_m_null_jmf.tab',out_file=null_ind,guessing=1005);
%import_files(input_file='!MCLAB/arbeitman_fru_network/exported_data_from_michelle/Repressed_Fru_m_null_jmf.tab',out_file=null_rep,guessing=713);


%macro remove_dups(in_data);

    proc sort data=&in_data(where=(Fbgn ^? ';')) nodupkey out=&in_data._nodups;
        by Fbgn;
        run;

    data &in_data._nodups;
        set &in_data._nodups;
        keep Fbgn;
        run;
%mend;

%remove_dups(null_ind); * 706 obs;
%remove_dups(null_rep); * 436 obs;

