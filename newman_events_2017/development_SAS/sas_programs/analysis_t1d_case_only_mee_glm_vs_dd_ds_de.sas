
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';


/* compare DD/DS genes with KW genes */

data kw_genes;
  set event.t1d_rank_glm_results_by_gene;
  keep gene_id flag_mee_p05;
run;

data dd_genes;
  set con.flag_genes_dd_ds_de;
  keep gene_id flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on flag_gene_on_multicell flag_gene_ds flag_gene_de;
run;

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

data multigene;
   set hg19.hg19_aceview_fusions_si_info;
   where flag_multigene=1;
   keep gene_id;
run;

proc sort data=kw_genes;
  by gene_id;
proc sort data=dd_genes;
  by gene_id;
proc sort data=ai nodup;
  by gene_id;
proc sort data=t1d nodup;
  by gene_id;
proc sort data=multigene nodup;
  by gene_id;
run;

data genes_to_compare;
  merge dd_genes (in=in1) kw_genes (in=in2) ai (in=in3) t1d (in=in4) multigene (in=in5);
  by gene_id;
  if in3 then flag_ai_gene=1; else flag_ai_gene=0;
  if in4 then flag_t1d_gene=1; else flag_t1d_gene=0;
  if in5 then flag_multigene=1; else flaG_multigene=0;
  if in1 and in2 then output;
run;

proc freq data=genes_to_compare noprint;
  tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*flag_gene_on_multicell*flag_gene_ds*flag_gene_de*flag_mee_p05 / out=gene_count;
run;

proc print data=gene_count;
run;

/*

                         flag_
flag_cd4_   flag_cd8_    cd19_     flag_gene_     flag_     flag_     flag_
 gene_on     gene_on    gene_on   on_multicell   gene_ds   gene_de   mee_p05   COUNT   PERCENT

    0           0          1            0           .         .         0          7     .
    1           1          1            1           0         1         0       2517   43.3816
    0           1          1            1           1         0         0          3    0.0517
    0           1          1            1           1         1         0          1    0.0172
    1           1          0            1           0         0         0          1    0.0172
    1           1          1            1           0         0         0        297    5.1189
    1           1          1            1           1         0         0          9    0.1551
    1           1          1            1           1         1         0       1310   22.5784
    0           1          1            1           0         0         0         13    0.2241
    0           0          1            0           .         .         1          1     .
    0           1          1            1           0         0         1          3    0.0517
    0           1          1            1           0         1         1          2    0.0345
    0           1          1            1           1         0         1         11    0.1896
    1           1          1            1           0         0         1          3    0.0517
    1           1          1            1           0         1         1        218    3.7573
    1           1          1            1           1         0         1          1    0.0172
    1           1          1            1           1         1         1       1413   24.3537

*/


proc freq data=genes_to_compare ;
  where flag_gene_on_multicell=1;
  tables flag_gene_ds*flag_mee_p05 ;
run;

/*
 Table of flag_gene_ds by flag_mee_p05

  flag_gene_ds     flag_mee_p05

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   2828 |    226 |   3054
           |  48.74 |   3.90 |  52.64
           |  92.60 |   7.40 |
           |  68.13 |  13.69 |
  ---------+--------+--------+
         1 |   1323 |   1425 |   2748
           |  22.80 |  24.56 |  47.36
           |  48.14 |  51.86 |
           |  31.87 |  86.31 |
  ---------+--------+--------+
  Total        4151     1651     5802
              71.54    28.46   100.00
*/

proc freq data=genes_to_compare ;
  where flag_gene_on_multicell=1;
  tables flag_gene_de*flag_mee_p05 ;
run;

/*
 Table of flag_gene_de by flag_mee_p05

  flag_gene_de     flag_mee_p05

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |    323 |     18 |    341
           |   5.57 |   0.31 |   5.88
           |  94.72 |   5.28 |
           |   7.78 |   1.09 |
  ---------+--------+--------+
         1 |   3828 |   1633 |   5461
           |  65.98 |  28.15 |  94.12
           |  70.10 |  29.90 |
           |  92.22 |  98.91 |
  ---------+--------+--------+
  Total        4151     1651     5802
              71.54    28.46   100.00
*/



proc freq data=genes_to_compare ;
  where flag_gene_on_multicell=1 and flag_ai_gene=1;
  tables flag_gene_ds*flag_mee_p05 ;
run;

/*
 flag_gene_ds     flag_mee_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |     66 |     26 |     92
          |  15.71 |   6.19 |  21.90
          |  71.74 |  28.26 |
          |  32.20 |  12.09 |
 ---------+--------+--------+
        1 |    139 |    189 |    328
          |  33.10 |  45.00 |  78.10
          |  42.38 |  57.62 |
          |  67.80 |  87.91 |
 ---------+--------+--------+
 Total         205      215      420
             48.81    51.19   100.00

*/

proc freq data=genes_to_compare ;
  where flag_gene_on_multicell=1  and flag_ai_gene=1;
  tables flag_gene_de*flag_mee_p05 ;
run;


/*

     flag_gene_de     flag_mee_p05

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            1 |    205 |    215 |    420
              |  48.81 |  51.19 | 100.00
              |  48.81 |  51.19 |
              | 100.00 | 100.00 |
     ---------+--------+--------+
     Total         205      215      420
                 48.81    51.19   100.00

*/


proc freq data=genes_to_compare ;
  where flag_gene_on_multicell=1 and flag_t1d_gene=1;
  tables flag_gene_ds*flag_mee_p05 ;
run;

/*
  flag_gene_ds     flag_mee_p05

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |     16 |      6 |     22
           |  16.33 |   6.12 |  22.45
           |  72.73 |  27.27 |
           |  29.09 |  13.95 |
  ---------+--------+--------+
         1 |     39 |     37 |     76
           |  39.80 |  37.76 |  77.55
           |  51.32 |  48.68 |
           |  70.91 |  86.05 |
  ---------+--------+--------+
  Total          55       43       98
              56.12    43.88   100.00


*/


proc freq data=genes_to_compare ;
  where flag_gene_on_multicell=1  and flag_t1d_gene=1;
  tables flag_gene_de*flag_mee_p05 ;
run;

/*

 flag_gene_de     flag_mee_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 |     55 |     43 |     98
          |  56.12 |  43.88 | 100.00
          |  56.12 |  43.88 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total          55       43       98
             56.12    43.88   100.00
*/


