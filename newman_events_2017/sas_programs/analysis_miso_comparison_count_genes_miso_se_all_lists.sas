ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Compare all gene lists for MISO */


data de_genes;
  set event.miso_refseq_exonskip_cmpr_de;
run;

data dtct_genes;
  set event.miso_refseq_exonskip_cmpr_dtct;
run;

data all_genes;
  set event.miso_refseq_exonskip_compare;
run;

proc sort data=all_genes;
   by gene_id;
proc sort data=dtct_genes;
   by gene_id;
proc sort data=de_genes;
   by gene_id;
run;

data all_miso2event_genes;
   merge all_genes dtct_genes de_genes;
   by gene_id;
run;

/* I want to count the genes with detected/analyzed SE events, but only for the set of Refseq genes */

proc freq data=all_miso2event_genes noprint;
   where flag_has_refseq=1;
   tables flag_has_miso_se_dtct*flag_exonskip_dtct_both*flag_exonskip_dtct_one / out=dtct_count;
run;

proc print data=dtct_count;
run;


/*
 flag_has_
  miso_se_    flag_exonskip_    flag_exonskip_
    dtct         dtct_both         dtct_one       COUNT

     1               0                 0            76
     1               0                 1           739
     1               1                 0             7
     1               1                 1           846


76 genes with a MISO SE event
739 genes with MISO SE and an ES junction in either NPC or OLD
7 genes with MISO SE and an ES junction in both NPC and OLD
846 genes with MISO, ES in both and ES in either
*/

/* I want to count the genes with DE/DD SE events, but only for the set of Refseq genes */

proc freq data=all_miso2event_genes;
   where flag_has_refseq=1 and flag_has_miso_se_dtct=1;
   tables flag_miso_se_diff_bf10*flag_event_se_dd
          flag_miso_se_diff_bf5*flag_event_se_dd   ;
run;



/*
 Table of flag_miso_se_diff_bf10 by flag_event_se_de

              flag_miso_se_diff_bf10
                        flag_event_se_de

              Frequency|
              Percent  |
              Row Pct  |
              Col Pct  |       0|  Total
              ---------+--------+
                     0 |   1656 |   1656
                       |  99.28 |  99.28
                       | 100.00 |
                       |  99.28 |
              ---------+--------+
                     1 |     12 |     12
                       |   0.72 |   0.72
                       | 100.00 |
                       |   0.72 |
              ---------+--------+
              Total        1668     1668
                         100.00   100.00

Table of flag_miso_se_diff_bf10 by flag_event_se_dd

        flag_miso_se_diff_bf10
                  flag_event_se_dd

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |   1039 |    617 |   1656
                 |  62.29 |  36.99 |  99.28
                 |  62.74 |  37.26 |
                 |  99.24 |  99.36 |
        ---------+--------+--------+
               1 |      8 |      4 |     12
                 |   0.48 |   0.24 |   0.72
                 |  66.67 |  33.33 |
                 |   0.76 |   0.64 |
        ---------+--------+--------+
        Total        1047      621     1668
                    62.77    37.23   100.00


 Table of flag_miso_se_diff_bf5 by flag_event_se_de

             flag_miso_se_diff_bf5
                       flag_event_se_de

             Frequency|
             Percent  |
             Row Pct  |
             Col Pct  |       0|  Total
             ---------+--------+
                    0 |   1650 |   1650
                      |  98.92 |  98.92
                      | 100.00 |
                      |  98.92 |
             ---------+--------+
                    1 |     18 |     18
                      |   1.08 |   1.08
                      | 100.00 |
                      |   1.08 |
             ---------+--------+
             Total        1668     1668
                        100.00   100.00

 Table of flag_miso_se_diff_bf5 by flag_event_se_dd

        flag_miso_se_diff_bf5
                  flag_event_se_dd

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |   1036 |    614 |   1650
                 |  62.11 |  36.81 |  98.92
                 |  62.79 |  37.21 |
                 |  98.95 |  98.87 |
        ---------+--------+--------+
               1 |     11 |      7 |     18
                 |   0.66 |   0.42 |   1.08
                 |  61.11 |  38.89 |
                 |   1.05 |   1.13 |
        ---------+--------+--------+
        Total        1047      621     1668
                    62.77    37.23   100.00
*/

proc freq data=all_miso2event_genes;
 where flag_has_refseq=1 and flag_has_miso_se_dtct=1 and (flag_exonskip_dtct_both=0 or flag_exonskip_dtct_one=0);
   tables           flag_miso_se_diff_bf10
                    flag_miso_se_diff_bf5;
run;

