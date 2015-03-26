/********************************************************************************
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';
* libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Import Motif Data */
    %macro import_motif(indata,outdata);

        proc import datafile=&indata out=&outdata dbms=CSV replace;
            getnames=yes;
            guessingrows=1000; 
            run;

        proc sort data=&outdata;
            by primary_fbgn;
            run;

    %mend import_motif;

    %import_motif('!MCLAB/arbeitman/arbeitman_Fru_network/motif_analysis/fru_a_results_up_and_down.csv',fru_a);
    %import_motif('!MCLAB/arbeitman/arbeitman_Fru_network/motif_analysis/fru_b_results_up_and_down.csv',fru_b);
    %import_motif('!MCLAB/arbeitman/arbeitman_Fru_network/motif_analysis/fru_c_results_up_and_down.csv',fru_c);

    %macro create_full_dataset(letter);

        data fru_&letter.2;
            set fru_&letter;
            rename motif_count = fru_&letter._motif_cnt;
            rename motif_positions = fru_&letter._motif_pos;
            rename motif_eval = fru_&letter._motif_eval;
            if motif_count = 0  then flag_fru_&letter._motif       = 0; else flag_fru_&letter._motif = 1;
            if motif_count = 0  then flag_fru_&letter._multi_motif = 0;
            if motif_count > 0  and  motif_count     <= 5  then flag_fru_&letter._multi_motif = 1;
            if motif_count > 5  and  motif_count     <= 10 then flag_fru_&letter._multi_motif = 2;
            if motif_count > 10 then flag_multi_motif = 3;
            keep primary_fbgn motif_count motif_positions motif_eval flag_fru_&letter._motif flag_fru_&letter._multi_motif;
            run;

        proc sort data=fru_&letter.2;
            by primary_fbgn;
            run;

    %mend create_full_dataset;

    %create_full_dataset(a); * 14903 obs;
    %create_full_dataset(b); * 14903 obs;
    %create_full_dataset(c); * 14903 obs;

    data FRU.motif_flags_and_cnts;
        merge fru_a2 (in=in1) fru_b2 (in=in2) fru_c2 (in=in3);
        by primary_fbgn;
        run;

/* Clean up */
    proc datasets nolist;
        delete fru_a fru_a2 fru_b fru_b2 fru_c fru_c2;
        run; quit;
