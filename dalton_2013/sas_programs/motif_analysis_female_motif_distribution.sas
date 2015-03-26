data flag_ind_rep;
    set FRU.flag_ind_rep;
    run;

data flag_motif;
    set FRU.motif_flags_and_cnts;
    keep primary_fbgn fru_a_motif_cnt fru_b_motif_cnt fru_c_motif_cnt;
    run;

proc sort data=flag_ind_rep;
    by primary_fbgn;
    run;

proc sort data=flag_motif;
    by primary_fbgn;
    run;

data ind_rep_motif;
    merge flag_ind_rep flag_motif;
    by primary_fbgn;
    keep primary_fbgn fru_a_motif_cnt fru_b_motif_cnt fru_c_motif_cnt 
    flag_female_a_ind 
    flag_female_b_ind 
    flag_female_c_ind 
    flag_female_a_rep 
    flag_female_b_rep 
    flag_female_c_rep; 
    run;

%macro out_loop(type);
    %macro letter_loop(letter);
        proc freq data=ind_rep_motif (where=(flag_female_&letter._&type = 1)) noprint; 
            table fru_&letter._motif_cnt /out=freq_&letter;
            run;

        data freq_&letter.2;
            set freq_&letter;
            rename fru_&letter._motif_cnt = motif_cnt;
            run;

        proc sort data=freq_&letter.2;
            by motif_cnt;
            run;

    %mend;
    %letter_loop(a);
    %letter_loop(b);
    %letter_loop(c);

    data freq_motif;
        set freq_a2 (in=in1) freq_b2 (in=in2) freq_c2 (in=in3) ;
        by motif_cnt;
        if in1 then motif = 'A';
        if in2 then motif = 'B';
        if in3 then motif = 'C';
        rename COUNT = num_genes;
        drop percent;
        run;

    proc export data=freq_motif
                outfile="!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/female_&type._motif_distribution.csv"
                dbms=CSV replace;
                putnames=yes;
                run;

%mend;
%out_loop(ind);
%out_loop(rep);

proc datasets nolist;
    delete FLAG_IND_REP FLAG_MOTIF FREQ_A FREQ_A2 FREQ_B FREQ_B2 FREQ_C FREQ_C2
    FREQ_MOTIF IND_REP_MOTIF;
    run;
