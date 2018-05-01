
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* Run rank test on complete T1D case only data (genes without multigene fusions only) */

data exon_rank;
  set eventloc.t1d_exons_mee_rank_cell_subj_v2;
  keep gene_id fusion_id subject_id rank_cd19 rank_cd8 rank_Cd4;
run;

proc transpose data=exon_rank out=exon_rank_stack;
   by gene_id fusion_id subject_id;
   var rank_cd19 rank_cd8 rank_cd4;
run;

data exon_rank_stack2;
   length cell_type $4.;
   set exon_rank_stack;
   if _NAME_ = "rank_cd19" then cell_type="CD19";
   if _NAME_ = "rank_cd8" then cell_type="CD8";
   if _NAME_ = "rank_cd4" then cell_type="CD4";
   rename COL1=exon_rank;
   drop _NAME_;
run;

ods listing close;
proc glm data=exon_rank_stack2 ;
  by gene_id;
  class cell_type fusion_id;
  model exon_rank = cell_type|fusion_id / ss3;
  output out=exon_rank_resid predicted=pred residual=resid student=stu;
  ods output ModelAnova=exon_rank_glm_type3;
run;
ods listing;

/* Transpose results and make permenant 
   I want to keep the DF, Type3 SS, Mean Seq, F-value and ProbF for each term */

data exon_rank_glm_type3_2;
  set exon_rank_glm_type3;
  length term $15.;
  if source="cell_type" then term="cell";
  else if source="fusion_id" then term="exon";
  else if source="cell_type*fusion_id" then term="cell_by_exon";
  keep gene_id term df ss ms fvalue probf;
run;

proc transpose data=exon_rank_glm_type3_2 out=df_sbys
     (rename=(cell=celltype_df exon=exon_df cell_by_exon=cell_by_exon_df) drop=_name_ _label_);
   by gene_id;
   var df;
   id term;
run;

proc transpose data=exon_rank_glm_type3_2 out=ss_sbys
     (rename=(cell=celltype_ss exon=exon_ss cell_by_exon=cell_by_exon_ss) drop=_name_ _label_);
   by gene_id;
   var ss;
   id term;
run;

proc transpose data=exon_rank_glm_type3_2 out=ms_sbys
     (rename=(cell=celltype_ms exon=exon_ms cell_by_exon=cell_by_exon_ms) drop=_name_ _label_);
   by gene_id;
   var ms;
   id term;
run;

proc transpose data=exon_rank_glm_type3_2 out=fvalue_sbys
     (rename=(cell=celltype_F exon=exon_F cell_by_exon=cell_by_exon_F) drop=_name_ _label_);
   by gene_id;
   var fvalue;
   id term;
run;


proc transpose data=exon_rank_glm_type3_2 out=probf_sbys
     (rename=(cell=celltype_P exon=exon_P cell_by_exon=cell_by_exon_P) drop=_name_ _label_);
   by gene_id;
   var probf;
   id term;
run;

data exon_rank_glm_sbys;
   merge df_sbys ss_sbys ms_sbys fvalue_sbys probf_sbys;
   by gene_id;
run;


/* Check: what is the distribution of significant diffs vs number of fusions? */

data flag_p05;
  set exon_rank_glm_sbys;
  if cell_by_exon_P ne . and cell_by_exon_P < 0.05 then flag_mee_p05=1; else flag_mee_p05=0;
run;


proc freq data=flag_p05;
  tables flag_mee_p05;
run;

/*
                                           Cumulative    Cumulative
  flag_mee_p05    Frequency     Percent     Frequency      Percent
  -----------------------------------------------------------------
             0        4158       71.57          4158        71.57
             1        1652       28.43          5810       100.00
*/

/* Make permenant */

data event.t1d_rank_glm_results_by_gene;
  set flag_p05;
run;


/* Check: are AI and T1D genes enriched? */

data ai t1d;
  set con.immunobase_gene_flags;
  if flag_immuno_gene=1 then output ai;
  if flag_immunobase_diabetes_gene=1 then output t1d;
run;

proc sort data=ai nodup;
  by gene_id;
proc sort data=t1d nodup;
  by gene_id;
proc sort data=flag_p05;
  by gene_id;
run;

data ai_t1d_check_p05;
  merge flag_p05 (in=in1) ai (in=in2) t1d (in=in3);
  by gene_id;
  if in2 then flag_ai_gene=1; else flag_ai_gene=0;
  if in3 then flag_t1d_gene=1; else flag_t1d_gene=0;
  if in1 then output;
run;

proc freq data=ai_t1d_check_p05;
  tables flag_mee_p05*flag_ai_gene / chisq;
run;

/*
  flag_mee_p05     flag_ai_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   3953 |    205 |   4158
           |  68.04 |   3.53 |  71.57
           |  95.07 |   4.93 |
           |  73.34 |  48.81 |
  ---------+--------+--------+
         1 |   1437 |    215 |   1652
           |  24.73 |   3.70 |  28.43
           |  86.99 |  13.01 |
           |  26.66 |  51.19 |
  ---------+--------+--------+
  Total        5390      420     5810
              92.77     7.23   100.00

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    115.2165    <.0001
   Likelihood Ratio Chi-Square    1    104.3614    <.0001
   Continuity Adj. Chi-Square     1    114.0141    <.0001
   Mantel-Haenszel Chi-Square     1    115.1966    <.0001
   Phi Coefficient                       0.1408
   Contingency Coefficient               0.1394
   Cramer's V                            0.1408


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      3953
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                     Sample Size = 5810
*/

proc freq data=ai_t1d_check_p05;
  tables flag_mee_p05*flag_t1d_gene / chisq;
run;

/*
 Table of flag_mee_p05 by flag_t1d_gene

  flag_mee_p05     flag_t1d_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   4103 |     55 |   4158
           |  70.62 |   0.95 |  71.57
           |  98.68 |   1.32 |
           |  71.83 |  56.12 |
  ---------+--------+--------+
         1 |   1609 |     43 |   1652
           |  27.69 |   0.74 |  28.43
           |  97.40 |   2.60 |
           |  28.17 |  43.88 |
  ---------+--------+--------+
  Total        5712       98     5810
              98.31     1.69   100.00

   Statistics for Table of flag_mee_p05 by flag_t1d_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     11.6837    0.0006
   Likelihood Ratio Chi-Square    1     10.7669    0.0010
   Continuity Adj. Chi-Square     1     10.9245    0.0009
   Mantel-Haenszel Chi-Square     1     11.6817    0.0006
   Phi Coefficient                       0.0448
   Contingency Coefficient               0.0448
   Cramer's V                            0.0448


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      4103
             Left-sided Pr <= F          0.9997
             Right-sided Pr >= F         0.0007

             Table Probability (P)       0.0004
             Two-sided Pr <= P           0.0010

                     Sample Size = 5810

*/


proc freq data=ai_t1d_check_p05;
  where flag_ai_gene=1;
  tables flag_mee_p05*flag_t1d_gene / chisq;
run;

/*
  Table of flag_mee_p05 by flag_t1d_gene

   flag_mee_p05     flag_t1d_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |    150 |     55 |    205
            |  35.71 |  13.10 |  48.81
            |  73.17 |  26.83 |
            |  46.58 |  56.12 |
   ---------+--------+--------+
          1 |    172 |     43 |    215
            |  40.95 |  10.24 |  51.19
            |  80.00 |  20.00 |
            |  53.42 |  43.88 |
   ---------+--------+--------+
   Total         322       98      420
               76.67    23.33   100.00

              The SAS System        16:19 Fr

    Statistics for Table of flag_mee_p05 by flag_t1d_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      2.7359    0.0981
    Likelihood Ratio Chi-Square    1      2.7392    0.0979
    Continuity Adj. Chi-Square     1      2.3675    0.1239
    Mantel-Haenszel Chi-Square     1      2.7294    0.0985
    Phi Coefficient                      -0.0807
    Contingency Coefficient               0.0804
    Cramer's V                           -0.0807

*/


