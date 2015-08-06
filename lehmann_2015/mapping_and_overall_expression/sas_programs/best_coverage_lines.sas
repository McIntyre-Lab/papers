proc import datafile='!MCLAB/cegs_sergey/pipeline_output/OE_alignment_summaries/aln_summary_fb551_non-redundant_fusions_20130912.csv' out=aln_summary dbms=csv replace;
    getnames=yes;
    run;

data aln_summary2;
    set aln_summary;
    aln_reads = bowtie1_aln + last_uniq;
    rename t_var_1 = line;
    rename t_var_2 = mating_status;
    rename t_var_3 = rep;
    run;

proc sort data=aln_summary2;
    by line;
    run;

proc means data=aln_summary2 noprint;
    by line;
    output out=sums sum(aln_reads)=sum_aln_reads;
    run;

proc sort data=sums;
    by DESCENDING sum_aln_reads;
    run;

data toplines;
    set sums;
    drop _type_ _freq_;
    run;

proc export data=toplines outfile='/home/jfear/mclab/cegs_sergey/reports/OE_alignments/all_lines_coverage.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Now keep Mating status separate */
proc sort data=aln_summary2;
    by line mating_status;
    run;

proc means data=aln_summary2 noprint;
    by line mating_status;
    output out=sums2 sum(aln_reads)=sum_aln_reads;
    run;

proc transpose data=sums2 out=flip;
    by line;
    id mating_status;
    var sum_aln_reads;
    run;

proc sort data=flip;
    by DESCENDING M V;
    run;

data toplines_MS;
    set flip;
    drop _name_;
    run;

proc export data=toplines_MS outfile='/home/jfear/mclab/cegs_sergey/reports/OE_alignments/all_lines_w_mating_status_coverage.csv' dbms=csv replace;
    putnames=yes;
    run;
