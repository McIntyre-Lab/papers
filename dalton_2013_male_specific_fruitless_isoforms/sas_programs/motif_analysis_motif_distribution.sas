proc freq data=FRU.motif_flags_and_cnts; 
    table fru_a_motif_cnt /out=freq_gene_w_motif_a;
    table fru_b_motif_cnt /out=freq_gene_w_motif_b;
    table fru_c_motif_cnt /out=freq_gene_w_motif_c;
    run;

proc sort data=freq_gene_w_motif_a;
    by motif_count;
    run;
proc sort data=freq_gene_w_motif_b;
    by motif_count;
    run;
proc sort data=freq_gene_w_motif_c;
    by motif_count;
    run;

data freq_motif;
    set freq_gene_w_motif_a (in=in1) freq_gene_w_motif_b (in=in2) freq_gene_w_motif_c (in=in3) ;
    by motif_count;
    if in1 then motif = 'A';
    if in2 then motif = 'B';
    if in3 then motif = 'C';
    rename COUNT = num_genes;
    drop percent;
    run;

proc export data=freq_motif
            outfile='!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_distribution.csv'
            dbms=CSV replace;
            putnames=yes;
            run;

