
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

/* I am testing 3 genes: IKZF3, FAM129C, IL8
   Using Wilcoxon rank, Krusal-Wallis and Cochran-Mantel-Haenszel tests
   to see which performs "best"/most appropriately for the data

   Expectation is that IKZF3 and FAM129C will be difference, IL8 will not be.

   For Wilcoxon rank and KW, analyze the pre-assigned ranks (1=MEE, 3=LEE, 2=others)

   For CHM, I think I need to count the number of times each exon is the MEE per group
   then run the CMH test using PROC FREQ. 
*/


/* Get my subset data */

data genes_for_test;
  length gene_id $12.;
  format gene_id $12.;
  informat gene_id $12.;
  set event.t1d_genes_for_mee_testing;
run;

/* Wilcoxon and KW tests -- transpose exon ranks and run with PROC NPAR1WAY
   For each I will print the result to ODS and paste it here
   */

data exon_rank;
  set genes_for_test;
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

proc sort data=exon_rank_stack2;
  by gene_id cell_type exon_rank;
run;

proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  class cell_type;
  var exon_rank;
  output out=cd4cd8cd19_kw_exon_rank WILCOXON;
run;

proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  where cell_type ne "CD8";
  class cell_type;
  var exon_rank;
  output out=cd4cd19_kw_exon_rank WILCOXON;
run;

proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  where cell_type ne "CD4";
  class cell_type;
  var exon_rank;
  output out=cd8cd19_kw_exon_rank WILCOXON;
run;

proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  where cell_type ne "CD19";
  class cell_type;
  var exon_rank;
  output out=cd4cd8_kw_exon_rank WILCOXON;
run;

proc print data=cd4cd8cd19_kw_exon_rank; quit;
proc print data=cd4cd8_kw_exon_rank; quit;
proc print data=cd4cd19_kw_exon_rank; quit;
proc print data=cd8cd19_kw_exon_rank; quit;



/* SUMMARY:

CD4 vs CD8 vs CD19:
Obs    gene_id         _VAR_         _KW_    DF_KW      P_KW
 1     FAM129C         exon_rank     279.41      2      0.00000
 2     FLT4            exon_rank    1970.42      2      0.00000
 3     IKZF3           exon_rank       0.07      2      0.96776
 4     IL8             exon_rank       0.00      2      1.00000


CD4 vs CD8:

gene_id  _VAR_          _WIL_   Z_WIL   PL_WIL      PR_WIL  P2_WIL   PTL_WIL     PTR_WIL  PT2_WIL     _KW_  DF_KW     P_KW
FAM129C  exon_rank  1709718.0  7.09569       .   6.4352E-13 0.00000        .  8.3234E-13  0.00000  50.3493      1  0.00000
FLT4     exon_rank  4289068.0  3.35684       .  0.000394199 0.00079        . 0.000397783  0.00080  11.2685      1  0.00079
IKZF3    exon_rank   505876.5  0.00000     0.5            . 1.00000      0.5          .   1.00000   0.0000      1  1.00000
IL8      exon_rank    25043.0  0.00000     0.5            . 1.00000      0.5          .   1.00000   0.0000      1  1.00000

CD4 vs CD19:

gene_id  _VAR_          _WIL_    Z_WIL   PL_WIL   PR_WIL   P2_WIL  PTL_WIL   PTR_WIL  PT2_WIL     _KW_  DF_KW      P_KW
FAM129C  exon_rank  1353033.0  -16.7779     0.0        .  0.00000      0.0         .  0.00000   281.50      1   0.00000
FLT4     exon_rank  5353593.0   34.6553       .  0.00000  0.00000        .   0.00000  0.00000  1200.99      1   0.00000
IKZF3    exon_rank   507140.5    0.2213       .  0.41245  0.82489        .   0.41246  0.82492     0.05      1   0.82482
IL8      exon_rank    25043.0    0.0000     0.5        .  1.00000      0.5         .  1.00000     0.00      1   1.00000

CD8 vs CD19:

gene_id  _VAR_          _WIL_    Z_WIL   PL_WIL   PR_WIL   P2_WIL  PTL_WIL   PTR_WIL  PT2_WIL     _KW_  DF_KW      P_KW
FAM129C  exon_rank  1464423.0  -10.1215     0.0        .  0.00000      0.0         .  0.00000   102.45      1   0.00000
FLT4     exon_rank  5422718.0   37.0691       .  0.00000  0.00000        .   0.00000  0.00000  1374.12      1   0.00000
IKZF3    exon_rank   507140.5    0.2213       .  0.41245  0.82489        .   0.41246  0.82492     0.05      1   0.82482
IL8      exon_rank    25043.0    0.0000     0.5        .  1.00000      0.5         .  1.00000     0.00      1   1.00000

From these, FAM129C and FLT4 are diff, IZKF3 and IL8 are not
*/

/* CMH test: prep data */

data exon_mee;
  set genes_for_test;
  keep gene_id fusion_id subject_id flag_mee_cd19 flag_mee_cd4 flag_mee_cd8;
run;

proc sort data=exon_mee;
   by gene_id fusion_id subject_id;
proc means data=exon_mee noprint;
   by gene_id fusion_id;
   var flag_mee_cd19 flag_mee_cd4 flag_mee_cd8;
   output out=mee_count_by_cell sum=;
run;

proc transpose data=mee_count_by_cell out=mee_count_stack;
   by gene_id fusion_id;
   var flag_mee_cd19 flag_mee_cd4 flag_mee_cd8;
run;

data mee_count_stack2;
  length cell_type $4.;
   set mee_count_stack;
   if _NAME_ = "flag_mee_cd19" then cell_type="CD19";
   if _NAME_ = "flag_mee_cd4" then cell_type="CD4";
   if _NAME_ = "flag_mee_cd8" then cell_type="CD8";
   drop _NAME_;
   rename col1=mee_count;
run;

proc freq data=mee_count_stack2;
   by gene_id;
   tables fusion_id*cell_type / cmh;
   weight mee_count;
run;

/* OUTPUT: skipping 3xN table output

FAM129C:

       Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

     Statistic    Alternative Hypothesis    DF       Value      Prob
     ---------------------------------------------------------------
         1        Nonzero Correlation        1     14.4125    0.0001
         2        Row Mean Scores Differ     2     28.2335    <.0001
         3        General Association       24     93.7943    <.0001
                    Total Sample Size = 237

FLT4:
            Summary Statistics for cell_type by fusion_id

     Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

   Statistic    Alternative Hypothesis    DF       Value      Prob
   ---------------------------------------------------------------
       1        Nonzero Correlation        1     10.0456    0.0015
       2        Row Mean Scores Differ     2     10.0714    0.0065
       3        General Association       44    108.5098    <.0001


                       Total Sample Size = 237


IKZF3:
         Summary Statistics for cell_type by fusion_id

  Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

Statistic    Alternative Hypothesis    DF       Value      Prob
---------------------------------------------------------------
    1        Nonzero Correlation        1      0.1424    0.7059
    2        Row Mean Scores Differ     2      1.3290    0.5145
    3        General Association        2      1.3290    0.5145


                    Total Sample Size = 237


IL8:
(no summary statistics -- MEE is the same for all cell types and subjects)

*/

/* CMH test: prep data */

data exon_rank;
  set genes_for_test;
  keep gene_id fusion_id subject_id rank_cd19 rank_cd4 rank_cd8;
run;

proc sort data=exon_rank;
   by gene_id fusion_id subject_id;
proc freq data=exon_rank noprint;
   by gene_id fusion_id;
   tables rank_cd4 / out=cd4_rank;
   tables rank_cd8 / out=cd8_rank;
   tables rank_cd19 / out=cd19_rank;
run;

data cd4_rank2;
   set cd4_rank;
   length cell_type $4.;
   cell_type="CD4";
   drop PERCENT;
   rename rank_cd4=exon_rank;
run;

data cd8_rank2;
   set cd8_rank;
   length cell_type $4.;
   cell_type="CD8";
   drop PERCENT;
   rename rank_cd8=exon_rank;
run;

data cd19_rank2;
   set cd19_rank;
   length cell_type $4.;
   cell_type="CD19";
   drop PERCENT;
   rename rank_cd19=exon_rank;
run;

data rank_stack;
   set cd4_rank2 cd8_rank2 cd19_rank2;
run;

proc sort data=rank_stack;
   by gene_id cell_type exon_rank fusion_id;
run;


proc freq data=rank_stack;
   by gene_id;
   tables cell_type*fusion_id*exon_rank / cmh  ;
   weight count;
   output out=cmh_stats cmh ;
run;



/* OUTPUT: skipping 3xN table output

FAM129C:


             Summary Statistics for fusion_id by exon_rank
                       Controlling for cell_type

      Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

    Statistic    Alternative Hypothesis    DF       Value      Prob
    ---------------------------------------------------------------
        1        Nonzero Correlation        1     34.3511    <.0001
        2        Row Mean Scores Differ    15    871.3420    <.0001
        3        General Association       30   1345.7082    <.0001


                        Total Sample Size = 3792

                             The SAS System         16:19 Friday, Octo


FLT4:

             Summary Statistics for fusion_id by exon_rank
                       Controlling for cell_type

      Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

    Statistic    Alternative Hypothesis    DF       Value      Prob
    ---------------------------------------------------------------
        1        Nonzero Correlation        1     82.9117    <.0001
        2        Row Mean Scores Differ    25    920.7026    <.0001
        3        General Association       50   1378.4741    <.0001


                        Total Sample Size = 6162
IKZF3:

           Summary Statistics for fusion_id by exon_rank
                     Controlling for cell_type

    Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

  Statistic    Alternative Hypothesis    DF       Value      Prob
  ---------------------------------------------------------------
      1        Nonzero Correlation        1    213.4438    <.0001
      2        Row Mean Scores Differ     8   1626.0020    <.0001
      3        General Association       16   3135.4155    <.0001


                      Total Sample Size = 2133

IL8:
           Summary Statistics for fusion_id by exon_rank
                     Controlling for cell_type

    Cochran-Mantel-Haenszel Statistics (Based on Table Scores)

  Statistic    Alternative Hypothesis    DF       Value      Prob
  ---------------------------------------------------------------
      1        Nonzero Correlation        1    471.0000    <.0001
      2        Row Mean Scores Differ     1    471.0000    <.0001
      3        General Association        1    471.0000    <.0001

  Total Sample Size = 474
*/

/*
http://www.lexjansen.com/wuss/2011/analy/Papers_Lai_G_74949.pdf
PROC CATMOD log-linear model as a k x R x C equivalent to the Breslow-Day k x 2 x 2 

data do not need to be collapsed!!

need gene_id, fusion_id, subject_id, cell_type, rank
*/

proc transpose data=exon_rank out=exon_rank_tall;
   by gene_id fusion_id subject_id;
   var rank_cd19 rank_cd8 rank_cd4;
run;

data exon_rank_tall2;
  set exon_rank_tall;
  length cell_type $4;
  if _NAME_ = "rank_cd19" then cell_type="CD19";
  else if _NAME_ = "rank_cd4" then cell_type="CD4";
  else if _NAME_ = "rank_cd8" then cell_type="CD8";
  else delete;
  drop _NAME_;
run;

proc sort data=exon_rank_tall2;
   by gene_id cell_type fusion_id exon_rank;
run;


proc catmod data=exon_rank_tall2;
   by gene_id;
   model cell_type * fusion_id * exon_rank = _response_ /
    ml noiter noresponse nodesign nogls NOPROFILE ZERO=sampling;
    loglin cell_type|fusion_id|exon_rank ;
run; quit;

