
/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* I want to now count the number of genes detected in each cell type pair that have at least 1 fusion or 1 splicing
   differentially detected. This will be evidence of alternative splice (gene must be detected in both cell types!)
*/

data splicing_summary;
   set con.gene_diff_splicing_summary;
run;

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=ai nodup;
  by gene_id;
proc sort data=t1d nodup;
  by gene_id;
proc sort data=splicing_summary;
  by gene_id;
run;

data splicing_summary_w_flags;
  merge splicing_summary (in=in1) ai (in=in2) t1d (in=in3);
  by gene_id;
  if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
  if in1 then output;
run;

proc freq data=splicing_summary_w_flags;
   tables flag_cd4cd8_event_dd flag_cd4cd19_event_dd flag_cd8cd19_event_dd flag_gene_any_dd;
run;

/*
     flag_cd4cd8_                             Cumulative    Cumulative
         event_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       36551       83.10         36551        83.10
                1        7434       16.90         43985       100.00

                        Frequency Missing = 3078


    flag_cd4cd19_                             Cumulative    Cumulative
         event_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33788       76.82         33788        76.82
                1       10197       23.18         43985       100.00

    flag_cd8cd19_                             Cumulative    Cumulative
         event_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33684       76.58         33684        76.58
                1       10301       23.42         43985       100.00

                        Frequency Missing = 3078


                                              Cumulative    Cumulative
 flag_gene_any_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       28961       65.84         28961        65.84
                1       15024       34.16         43985       100.00

*/

proc freq data=splicing_summary_w_flags;
   where flag_immuno_gene=1;
   tables flag_cd4cd8_event_dd flag_cd4cd19_event_dd flag_cd8cd19_event_dd flag_gene_any_dd;
run;

/*

      flag_cd4cd8_                             Cumulative    Cumulative
          event_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0         867       52.20           867        52.20
                 1         794       47.80          1661       100.00

                          Frequency Missing = 29


     flag_cd4cd19_                             Cumulative    Cumulative
          event_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0         649       39.07           649        39.07
                 1        1012       60.93          1661       100.00

     flag_cd8cd19_                             Cumulative    Cumulative
          event_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0         642       38.65           642        38.65
                 1        1019       61.35          1661       100.00

                          Frequency Missing = 29


                                               Cumulative    Cumulative
  flag_gene_any_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0         381       22.94           381        22.94
                 1        1280       77.06          1661       100.00

*/

proc freq data=splicing_summary_w_flags;
   where flag_immunobase_diabetes_gene=1;
   tables flag_cd4cd8_event_dd flag_cd4cd19_event_dd flag_cd8cd19_event_dd flag_gene_any_dd;
run;

/*
       flag_cd4cd8_                             Cumulative    Cumulative
           event_dd    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0         215       54.02           215        54.02
                  1         183       45.98           398       100.00

                           Frequency Missing = 7


      flag_cd4cd19_                             Cumulative    Cumulative
           event_dd    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0         160       40.20           160        40.20
                  1         238       59.80           398       100.00

       flag_cd8cd19_                             Cumulative    Cumulative
            event_dd    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0         158       39.70           158        39.70
                   1         240       60.30           398       100.00

                            Frequency Missing = 7


                                                 Cumulative    Cumulative
    flag_gene_any_dd    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0         105       26.38           105        26.38
                   1         293       73.62           398       100.00
*/


