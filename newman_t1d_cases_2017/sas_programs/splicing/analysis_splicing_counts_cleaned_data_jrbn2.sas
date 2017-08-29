/* Redo counts and check proportions of DE and DD splicing events - cleaned dataset

Here I am going to look at the three DE/DD scenarios again:

1. DE events expressed in all three tissues
2. DE events expressed in all three tissues PLUS tissue-specific events
3. DE events expressed in all three tissues PLUS two-tissue events

*/

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';


/* Using cleaned dataset */

data immunogenes;
   set con.immunogene_flags;
run;

proc sort data=immunogenes;
   by gene_id;
proc sort data=splicing.splicing_results_clean;
   by gene_id;
run;

data splicing_results_clean_w_flags;
  merge splicing.splicing_results_clean (in=in1) immunogenes (in=in2);
   by gene_id;
  if in1 and in2 then output;
  else if in1 then do;
      flag_pseudogene=0;
      flag_immuno_gene=0;
      flag_diabetes_gene=0;
  end;
run;

data il21_splicing;
   set splicing.splicing_results_clean;
   if gene_id='IL21';
run;



/* Get total counts */

proc freq data=splicing_results_clean_w_flags;
   tables flag_immuno_gene*flag_diabetes_gene;
run;

* 5856526 possible (genes with expression);
* 645870 possible (autoimmune genes with expression);
* 8177 possible (T1D genes with expression);

/* Get counts - all genes */

* total splicing events;
proc freq data=splicing_results_clean_w_flags noprint;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_total;
run;

* exon skipping events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_exonskip=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_es;
run;

* alt donor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_ad;
run;

* alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_aa;
run;

* alt donor and alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_ada;
run;

* intron retention events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_intron_retention=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_ir;
run;

* annotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_annot;
run;

* unannotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=all_splicing_unannot;
run;

/*
	OFF	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	CD4/8/19	TOTAL	SINGLE	DOUBLE	%SINGLE	%COMMON
TOTAL	5646639	4349	6478	10008	15961	1913	2973	168205		209887	20835	20847	9.93%	80.14%
ES	5081166	1241	1732	2460	4267	471	721	35189		46081	5433	5459	11.79%	76.36%
AD	1088974	430	515	893	1480	182	245	17407		21152	1838	1907	8.69%	82.29%
AA	1101190	403	521	943	1475	203	261	18070		21876	1867	1939	8.53%	82.60%
ADA	788805	554	525	1016	1952	209	319	24131		28706	2095	2480	7.30%	84.06%
IR	184288	1417	2432	2778	4777	509	1058	27921		40892	6627	6344	16.21%	68.28%
Annot	131703	1852	2540	5182	7460	997	1313	110934		130278	9574	9770	7.35%	85.15%
Unannot	5330648	1080	1506	2048	3724	407	602	29350		38717	4634	4733	11.97%	75.81%


*/

/* Get counts - autoimmune genes */

* total splicing events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_total;
run;

* exon skipping events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_es;
run;

* alt donor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_ad;
run;

* alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_aa;
run;

* alt donor and alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_ada;
run;

* intron retention events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_ir;
run;

* annotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=1 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_annot;
run;

* unannotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=immuno_splicing_unannot;
run;


/*
	OFF	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	CD4/8/19	TOTAL	SINGLE	DOUBLE	%SINGLE	%COMMON
TOTAL	621226	483	675	1232	1935	219	334	19766		24644	2390	2488	9.70%	80.21%
ES	566184	157	141	310	539	53	73	3689		4962	608	665	12.25%	74.35%
AD	126090	43	50	125	197	17	27	2072		2531	218	241	8.61%	81.86%
AA	128323	51	52	127	173	24	29	2179		2635	230	226	8.73%	82.69%
ADA	109898	71	56	162	304	26	43	3599		4261	289	373	6.78%	84.46%
IR	16180	170	290	386	584	68	131	3065		4694	846	783	18.02%	65.30%
Annot	10625	176	270	580	894	95	143	13337		15495	1026	1132	6.62%	86.07%
Unannot	594421	137	115	266	457	56	60	3364		4455	518	573	11.63%	75.51%



*/


/* Get counts - diabetes genes */

* total splicing events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_total;
run;

* exon skipping events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_exonskip=1 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_es;
run;

* alt donor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_ad;
run;

* alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_aa;
run;

* alt donor and alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_ada;
run;

* intron retention events;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_intron_retention=1 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_ir;
run;

* annotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=1 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_annot;
run;

* unannotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=t1d_splicing_unannot;
run;


/*
	OFF	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	CD4/8/19
TOTAL	7713	9	22	63	61	11	5	293
ES	6804	3	5	18	19	1	0	65
AD	1536	2	1	9	5	1	0	27
AA	1731	0	3	10	3	3	0	34
ADA	1190	1	0	7	6	0	0	34
IR	261	4	14	15	24	3	5	47
Annot	187	2	2	35	15	8	0	204
Unannot	7265	3	6	13	22	0	0	42


*/

/* Get DE and DD counts */

data splicing_clean_w_de_flags;
  set splicing_results_clean_w_flags;

  if flag_CD19_on=1 or flag_CD4_on=1 or flag_CD8_on=1 then flag_any_on=1;
  else flag_any_on=0;

  if flag_CD19_on=1 and flag_CD4_on=1 and flag_CD8_on=1 then flag_all_on=1;
  else flag_all_on=0;

  /* Tissue specificity */

  if flag_CD19_on=1 and flag_CD4_on=0 and flag_CD8_on=0 then flag_tissue_specific=1;
  else if flag_CD19_on=0 and flag_CD4_on=1 and flag_CD8_on=0 then flag_tissue_specific=1;
  else if flag_CD19_on=0 and flag_CD4_on=0 and flag_CD8_on=1 then flag_tissue_specificn=1;
  else flag_tissue_specific=0;

  /* Scenario 1 - only DE events expressed in all 3 tissues */
  
  if flag_all_on=1 and flag_anova_fdr_05=1 then flag_event_de=1;
  else if flag_all_on=1 and flag_anova_fdr_05=0 then flag_event_de=0;
  else if flag_all_on=0 then flag_event_de=.;

  /* Scenario 2 - only DE events expressed in all 3 tissues PLUS tissue-specific events */

  if flag_all_on=1 then do;
     if flag_anova_fdr_05=1 then flag_event_dd_1=1;
     else flag_event_dd_1=0;
     end;
  else if flag_all_on=0 then do;
      if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_event_dd_1=1;
      else if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=0 then flag_event_dd_1=1;
      else if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=1 then flag_event_dd_1=1;
      else flag_event_dd_1=.;
      end;

  /* Scenario 3 - only DE events expressed in all 3 tissues PLUS two-tissue events */


  if flag_all_on=1 then do;
     if flag_anova_fdr_05=1 then flag_event_dd_2=1;
     else flag_event_dd_2=0;
     end;
  else if flag_all_on=0 then do;
      if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_event_dd_2=1;
      else flag_event_dd_2=.;
      end;

run;

/* Get DE and DD counts - all genes */


* total splicing events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* exon skipping events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_exonskip=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor and alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* intron retention events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* annotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* unannotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


/*
	Scenario 1	Scenario 2	Scenario 3
	NOT_DE	DE	NOT_DD1	DD1	NOT_DD2	DD2
TOTAL	76835	91370	76835	112205	76835	133052
ES	16410	18779	16410	24212	16410	29671
AD	7332	10075	7332	11913	7332	13820
AA	7722	10348	7722	12215	7722	14154
ADA	9949	14182	9949	16277	9949	18757
IR	17256	10665	17256	17292	17256	23636
Annot	45820	65114	45820	74688	45820	84458
Unannot	13759	15591	13759	20225	13759	24958

*/


/* Get DE and DD counts - autoimmune genes */


* total splicing events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* exon skipping events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor and alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* intron retention events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* annotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=1 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* unannotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

/*
	Scenario 1	Scenario 2	Scenario 3
	NOT_DE	DE	NOT_DD1	DD1	NOT_DD2	DD2
TOTAL	8343	11423	8343	13813	8343	16301
ES	1553	2136	1553	2744	1553	3409
AD	815	1257	815	1475	815	1716
AA	814	1365	814	1595	814	1821
ADA	1378	2221	1378	2510	1378	2883
IR	1906	1159	1906	2005	1906	2788
Annot	5037	8300	5037	9326	5037	10458
Unannot	1400	1964	1400	2482	1400	3055
*/



/* Get DE and DD counts - diabetes genes */


* total splicing events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* exon skipping events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_exonskip=1 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor and alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* intron retention events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=1 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* annotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=1 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* unannotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=0 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


/*
	Scenario 1	Scenario 2	Scenario 3
	NOT_DE	DE	NOT_DD1	DD1	NOT_DD2	DD2
TOTAL	88	205	88	299	88	376
ES	22	43	22	69	22	89
AD	7	20	7	32	7	38
AA	6	28	6	41	6	47
ADA	14	20	14	28	14	34
IR	23	24	23	57	23	89
Annot	52	152	52	191	52	214
Unannot	13	29	13	51	13	73

*/

/* Now need to look at proportions:

autoimmune among all detected
diabetes among autoimmune detected

autoimmune among de/dd detected
diabetes among de/dd autoimmune detected

*/

/****** OVERALL *******/

/* Autoimmune among all detected */

proc freq data=splicing_clean_w_de_flags;
   tables flag_any_on*flag_immuno_gene / chisq ;
run;

/*
 Table of flag_any_on by flag_immuno_gene

   flag_any_on     flag_immuno_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |5025413 | 621226 |5646639
            |  85.81 |  10.61 |  96.42
            |  89.00 |  11.00 |
            |  96.44 |  96.18 |
   ---------+--------+--------+
          1 | 185243 |  24644 | 209887
            |   3.16 |   0.42 |   3.58
            |  88.26 |  11.74 |
            |   3.56 |   3.82 |
   ---------+--------+--------+
   Total     5210656   645870  5856526
               88.97    11.03   100.00

 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1    112.8960    <.0001
  Likelihood Ratio Chi-Square    1    110.9080    <.0001
  Continuity Adj. Chi-Square     1    112.8206    <.0001
  Mantel-Haenszel Chi-Square     1    112.8960    <.0001
  Phi Coefficient                       0.0044
  Contingency Coefficient               0.0044
  Cramer's V                            0.0044


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    5025413
           Left-sided Pr <= F           1.0000
           Right-sided Pr >= F          <.0001

           Table Probability (P)        <.0001
           Two-sided Pr <= P            <.0001

                  Sample Size = 5856526

*/


/* T1D among autoimmune detected */

proc freq data=splicing_clean_w_de_flags;
   where flag_immuno_gene=1;
   tables flag_any_on*flag_diabetes_gene / chisq ;
run;

/*
Table of flag_any_on by flag_diabetes_gene

   flag_any_on     flag_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 | 613513 |   7713 | 621226
            |  94.99 |   1.19 |  96.18
            |  98.76 |   1.24 |
            |  96.21 |  94.33 |
   ---------+--------+--------+
          1 |  24180 |    464 |  24644
            |   3.74 |   0.07 |   3.82
            |  98.12 |   1.88 |
            |   3.79 |   5.67 |
   ---------+--------+--------+
   Total      637693     8177   645870
               98.73     1.27   100.00


Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     77.9711    <.0001
  Likelihood Ratio Chi-Square    1     68.2477    <.0001
  Continuity Adj. Chi-Square     1     77.4589    <.0001
  Mantel-Haenszel Chi-Square     1     77.9709    <.0001
  Phi Coefficient                       0.0110
  Contingency Coefficient               0.0110
  Cramer's V                            0.0110


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    613513
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 645870

*/


/* Autoimmune among all DE - Scenario 1 */

proc freq data=splicing_clean_w_de_flags;
   where flag_all_on=1;
   tables flag_event_de*flag_immuno_gene / chisq ;
run;

/*
  Table of flag_event_de by flag_immuno_gen

    flag_event_de     flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  68492 |   8343 |  76835
             |  40.72 |   4.96 |  45.68
             |  89.14 |  10.86 |
             |  46.14 |  42.21 |
    ---------+--------+--------+
           1 |  79947 |  11423 |  91370
             |  47.53 |   6.79 |  54.32
             |  87.50 |  12.50 |
             |  53.86 |  57.79 |
    ---------+--------+--------+
    Total      148439    19766   168205
                88.25    11.75   100.00

               The SAS System         11:4

 Statistics for Table of flag_event_de by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    108.7221    <.0001
   Likelihood Ratio Chi-Square    1    109.1832    <.0001
   Continuity Adj. Chi-Square     1    108.5637    <.0001
   Mantel-Haenszel Chi-Square     1    108.7215    <.0001
   Phi Coefficient                       0.0254
   Contingency Coefficient               0.0254
   Cramer's V                            0.0254


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     68492
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 168205


*/


/* Autoimmune among all DE - Scenario 2 */

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_event_dd_1*flag_immuno_gene / chisq ;
run;

/*
Table of flag_event_dd_1 by flag_immuno_gene

    flag_event_dd_1
              flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  68492 |   8343 |  76835
             |  36.23 |   4.41 |  40.64
             |  89.14 |  10.86 |
             |  41.04 |  37.66 |
    ---------+--------+--------+
           1 |  98392 |  13813 | 112205
             |  52.05 |   7.31 |  59.36
             |  87.69 |  12.31 |
             |  58.96 |  62.34 |
    ---------+--------+--------+
    Total      166884    22156   189040
                88.28    11.72   100.00


Statistics for Table of flag_event_dd_1 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     92.9510    <.0001
   Likelihood Ratio Chi-Square    1     93.6337    <.0001
   Continuity Adj. Chi-Square     1     92.8107    <.0001
   Mantel-Haenszel Chi-Square     1     92.9505    <.0001
   Phi Coefficient                       0.0222
   Contingency Coefficient               0.0222
   Cramer's V                            0.0222


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     68492
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

               Effective Sample Size = 189040
                 Frequency Missing = 20847
*/


/* Autoimmune among all DE - Scenario 3 */

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq ;
run;

/*
   Table of flag_event_dd_2 by flag_immuno_gene

       flag_event_dd_2
                 flag_immuno_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |  68492 |   8343 |  76835
                |  32.63 |   3.97 |  36.61
                |  89.14 |  10.86 |
                |  36.97 |  33.85 |
       ---------+--------+--------+
              1 | 116751 |  16301 | 133052
                |  55.63 |   7.77 |  63.39
                |  87.75 |  12.25 |
                |  63.03 |  66.15 |
       ---------+--------+--------+
       Total      185243    24644   209887
                   88.26    11.74   100.00

                  The SAS System         11:48 Mon

Statistics for Table of flag_event_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     91.2394    <.0001
   Likelihood Ratio Chi-Square    1     92.1370    <.0001
   Continuity Adj. Chi-Square     1     91.1050    <.0001
   Mantel-Haenszel Chi-Square     1     91.2390    <.0001
   Phi Coefficient                       0.0208
   Contingency Coefficient               0.0208
   Cramer's V                            0.0208


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     68492
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 209887


*/



/* Diabetes among autoimmune DE - Scenario 1 */

proc freq data=splicing_clean_w_de_flags;
   where flag_all_on=1 and flag_immuno_gene=1;
   tables flag_event_de*flag_diabetes_gene / chisq ;
run;

/*
  Table of flag_event_de by flag_diabetes_gene

      flag_event_de
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   8255 |     88 |   8343
               |  41.76 |   0.45 |  42.21
               |  98.95 |   1.05 |
               |  42.39 |  30.03 |
      ---------+--------+--------+
             1 |  11218 |    205 |  11423
               |  56.75 |   1.04 |  57.79
               |  98.21 |   1.79 |
               |  57.61 |  69.97 |
      ---------+--------+--------+
      Total       19473      293    19766
                  98.52     1.48   100.00

 Statistics for Table of flag_event_de by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     18.0720    <.0001
   Likelihood Ratio Chi-Square    1     18.7562    <.0001
   Continuity Adj. Chi-Square     1     17.5689    <.0001
   Mantel-Haenszel Chi-Square     1     18.0711    <.0001
   Phi Coefficient                       0.0302
   Contingency Coefficient               0.0302
   Cramer's V                            0.0302


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      8255
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 19766
*/


/* Diabetes among autoimmune DE - Scenario 2 */

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_1*flag_diabetes_gene / chisq ;
run;

/*

Table of flag_event_dd_1 by flag_diabetes_gene

     flag_event_dd_1
               flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   8255 |     88 |   8343
              |  37.26 |   0.40 |  37.66
              |  98.95 |   1.05 |
              |  37.92 |  22.74 |
     ---------+--------+--------+
            1 |  13514 |    299 |  13813
              |  60.99 |   1.35 |  62.34
              |  97.84 |   2.16 |
              |  62.08 |  77.26 |
     ---------+--------+--------+
     Total       21769      387    22156
                 98.25     1.75   100.00

 Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     37.3320    <.0001
     Likelihood Ratio Chi-Square    1     40.1578    <.0001
     Continuity Adj. Chi-Square     1     36.6881    <.0001
     Mantel-Haenszel Chi-Square     1     37.3303    <.0001
     Phi Coefficient                       0.0410
     Contingency Coefficient               0.0410
     Cramer's V                            0.0410


 Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      8255
               Left-sided Pr <= F          1.0000
               Right-sided Pr >= F         <.0001

               Table Probability (P)       <.0001
               Two-sided Pr <= P           <.0001

                 Effective Sample Size = 22156

*/


/* Diabetes among autoimmune DE - Scenario 3 */

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq ;
run;

/*
Table of flag_event_dd_2 by flag_diabetes_gene

     flag_event_dd_2
               flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   8255 |     88 |   8343
              |  33.50 |   0.36 |  33.85
              |  98.95 |   1.05 |
              |  34.14 |  18.97 |
     ---------+--------+--------+
            1 |  15925 |    376 |  16301
              |  64.62 |   1.53 |  66.15
              |  97.69 |   2.31 |
              |  65.86 |  81.03 |
     ---------+--------+--------+
     Total       24180      464    24644
                 98.12     1.88   100.00

                The SAS System         11:48 Mon


Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     46.8128    <.0001
    Likelihood Ratio Chi-Square    1     51.5640    <.0001
    Continuity Adj. Chi-Square     1     46.1376    <.0001
    Mantel-Haenszel Chi-Square     1     46.8109    <.0001
    Phi Coefficient                       0.0436
    Contingency Coefficient               0.0435
    Cramer's V                            0.0436


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      8255
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 24644

*/




/* Proportion of detected ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1;
   tables flag_any_on*flag_immuno_gene / chisq;
run;


/*
Table of flag_any_on by flag_immuno_gene

  flag_any_on     flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |4514982 | 566184 |5081166
           |  88.06 |  11.04 |  99.10
           |  88.86 |  11.14 |
           |  99.10 |  99.13 |
  ---------+--------+--------+
         1 |  41119 |   4962 |  46081
           |   0.80 |   0.10 |   0.90
           |  89.23 |  10.77 |
           |   0.90 |   0.87 |
  ---------+--------+--------+
  Total     4556101   571146  5127247
              88.86    11.14   100.00


 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      6.4808    0.0109
  Likelihood Ratio Chi-Square    1      6.5443    0.0105
  Continuity Adj. Chi-Square     1      6.4430    0.0111
  Mantel-Haenszel Chi-Square     1      6.4808    0.0109
  Phi Coefficient                      -0.0011
  Contingency Coefficient               0.0011
  Cramer's V                           -0.0011


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    4514982
           Left-sided Pr <= F           0.0054
           Right-sided Pr >= F          0.9948

           Table Probability (P)        0.0002
           Two-sided Pr <= P            0.0107

                  Sample Size = 5127247





*/

/* Proportion of DE ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_all_on=1;
   tables flag_event_de*flag_immuno_gene / chisq;
run;

/*
Table of flag_event_de by flag_immuno_gene

   flag_event_de     flag_immuno_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  14857 |   1553 |  16410
            |  42.22 |   4.41 |  46.63
            |  90.54 |   9.46 |
            |  47.17 |  42.10 |
   ---------+--------+--------+
          1 |  16643 |   2136 |  18779
            |  47.30 |   6.07 |  53.37
            |  88.63 |  11.37 |
            |  52.83 |  57.90 |
   ---------+--------+--------+
   Total       31500     3689    35189
               89.52    10.48   100.00

Statistics for Table of flag_event_de by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     34.0673    <.0001
  Likelihood Ratio Chi-Square    1     34.2340    <.0001
  Continuity Adj. Chi-Square     1     33.8640    <.0001
  Mantel-Haenszel Chi-Square     1     34.0664    <.0001
  Phi Coefficient                       0.0311
  Contingency Coefficient               0.0311
  Cramer's V                            0.0311


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     14857
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 35189
*/


/* Proportion of DD_1 ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1;
   tables flag_event_dd_1*flag_immuno_gene / chisq;
run;

/*
Table of flag_event_dd_1 by flag_immuno_gene

    flag_event_dd_1
              flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  14857 |   1553 |  16410
             |  36.57 |   3.82 |  40.40
             |  90.54 |   9.46 |
             |  40.90 |  36.14 |
    ---------+--------+--------+
           1 |  21468 |   2744 |  24212
             |  52.85 |   6.75 |  59.60
             |  88.67 |  11.33 |
             |  59.10 |  63.86 |
    ---------+--------+--------+
    Total       36325     4297    40622
                89.42    10.58   100.00

Statistics for Table of flag_event_dd_1 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     36.1387    <.0001
   Likelihood Ratio Chi-Square    1     36.5479    <.0001
   Continuity Adj. Chi-Square     1     35.9413    <.0001
   Mantel-Haenszel Chi-Square     1     36.1378    <.0001
   Phi Coefficient                       0.0298
   Contingency Coefficient               0.0298
   Cramer's V                            0.0298

         Fisher's Exact Test
  ----------------------------------
  Cell (1,1) Frequency (F)     14857
  Left-sided Pr <= F          1.0000
  Right-sided Pr >= F         <.0001

  Table Probability (P)       <.0001
  Two-sided Pr <= P           <.0001

    Effective Sample Size = 40622
*/


/* Proportion of DD_2 ES events that are autoimmune events */


proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;


/*
  Table of flag_event_dd_2 by flag_immuno_gene

      flag_event_dd_2
                flag_immuno_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  14857 |   1553 |  16410
               |  32.24 |   3.37 |  35.61
               |  90.54 |   9.46 |
               |  36.13 |  31.30 |
      ---------+--------+--------+
             1 |  26262 |   3409 |  29671
               |  56.99 |   7.40 |  64.39
               |  88.51 |  11.49 |
               |  63.87 |  68.70 |
      ---------+--------+--------+
      Total       41119     4962    46081
                  89.23    10.77   100.00


Statistics for Table of flag_event_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     45.1198    <.0001
   Likelihood Ratio Chi-Square    1     45.9146    <.0001
   Continuity Adj. Chi-Square     1     44.9092    <.0001
   Mantel-Haenszel Chi-Square     1     45.1188    <.0001
   Phi Coefficient                       0.0313
   Contingency Coefficient               0.0313
   Cramer's V                            0.0313


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     14857
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46081

*/





/* Proportion of detected IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1;
   tables flag_any_on*flag_immuno_gene / chisq;
run;


/*

Table of flag_any_on by flag_immuno_gene

  flag_any_on     flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 168108 |  16180 | 184288
           |  74.65 |   7.19 |  81.84
           |  91.22 |   8.78 |
           |  82.28 |  77.51 |
  ---------+--------+--------+
         1 |  36198 |   4694 |  40892
           |  16.08 |   2.08 |  18.16
           |  88.52 |  11.48 |
           |  17.72 |  22.49 |
  ---------+--------+--------+
  Total      204306    20874   225180
              90.73     9.27   100.00

Statistics for Table of flag_any_on by flag_immuno_gene

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1    289.9178    <.0001
 Likelihood Ratio Chi-Square    1    275.6400    <.0001
 Continuity Adj. Chi-Square     1    289.5970    <.0001
 Mantel-Haenszel Chi-Square     1    289.9165    <.0001
 Phi Coefficient                       0.0359
 Contingency Coefficient               0.0359
 Cramer's V                            0.0359


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)    168108
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 225180


*/

/* Proportion of DE IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_all_on=1;
   tables flag_event_de*flag_immuno_gene / chisq;
run;

/*
   Table of flag_event_de by flag_immuno_gene

      flag_event_de     flag_immuno_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  15350 |   1906 |  17256
               |  54.98 |   6.83 |  61.80
               |  88.95 |  11.05 |
               |  61.76 |  62.19 |
      ---------+--------+--------+
             1 |   9506 |   1159 |  10665
               |  34.05 |   4.15 |  38.20
               |  89.13 |  10.87 |
               |  38.24 |  37.81 |
      ---------+--------+--------+
      Total       24856     3065    27921
                  89.02    10.98   100.00

Statistics for Table of flag_event_de by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      0.2140    0.6437
  Likelihood Ratio Chi-Square    1      0.2142    0.6435
  Continuity Adj. Chi-Square     1      0.1961    0.6579
  Mantel-Haenszel Chi-Square     1      0.2140    0.6437
  Phi Coefficient                      -0.0028
  Contingency Coefficient               0.0028
  Cramer's V                           -0.0028


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     15350
            Left-sided Pr <= F          0.3293
            Right-sided Pr >= F         0.6849

            Table Probability (P)       0.0141
            Two-sided Pr <= P           0.6505

                   Sample Size = 27921


*/


/* Proportion of DD_1 IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1;
   tables flag_event_dd_1*flag_immuno_gene / chisq;
run;

/*
   Table of flag_event_dd_1 by flag_immuno_gene

       flag_event_dd_1
                 flag_immuno_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |  15350 |   1906 |  17256
                |  44.43 |   5.52 |  49.95
                |  88.95 |  11.05 |
                |  50.10 |  48.73 |
       ---------+--------+--------+
              1 |  15287 |   2005 |  17292
                |  44.25 |   5.80 |  50.05
                |  88.41 |  11.59 |
                |  49.90 |  51.27 |
       ---------+--------+--------+
       Total       30637     3911    34548
                   88.68    11.32   100.00

  Statistics for Table of flag_event_dd_1 by flag_immuno_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      2.5980    0.1070
     Likelihood Ratio Chi-Square    1      2.5983    0.1070
     Continuity Adj. Chi-Square     1      2.5436    0.1107
     Mantel-Haenszel Chi-Square     1      2.5980    0.1070
     Phi Coefficient                       0.0087
     Contingency Coefficient               0.0087
     Cramer's V                            0.0087

 Statistics for Table of flag_event_dd_1 by flag_immuno_gene

                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     15350
              Left-sided Pr <= F          0.9483
              Right-sided Pr >= F         0.0554

              Table Probability (P)       0.0037
              Two-sided Pr <= P           0.1104

                Effective Sample Size = 34548

*/


/* Proportion of DD_2 IR events that are autoimmune events */


proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;


/*
  Table of flag_event_dd_2 by flag_immuno_gene

      flag_event_dd_2
                flag_immuno_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  15350 |   1906 |  17256
               |  37.54 |   4.66 |  42.20
               |  88.95 |  11.05 |
               |  42.41 |  40.61 |
      ---------+--------+--------+
             1 |  20848 |   2788 |  23636
               |  50.98 |   6.82 |  57.80
               |  88.20 |  11.80 |
               |  57.59 |  59.39 |
      ---------+--------+--------+
      Total       36198     4694    40892
                  88.52    11.48   100.00

Statistics for Table of flag_event_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      5.5233    0.0188
   Likelihood Ratio Chi-Square    1      5.5406    0.0186
   Continuity Adj. Chi-Square     1      5.4498    0.0196
   Mantel-Haenszel Chi-Square     1      5.5232    0.0188
   Phi Coefficient                       0.0116
   Contingency Coefficient               0.0116
   Cramer's V                            0.0116


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     15350
             Left-sided Pr <= F          0.9911
             Right-sided Pr >= F         0.0097

             Table Probability (P)       0.0008
             Two-sided Pr <= P           0.0193

                    Sample Size = 40892

*/





/* Proportion of detected autoimmune ES events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_any_on*flag_diabetes_gene / chisq;
run;


/*
     Table of flag_any_on by flag_diabetes_gene

        flag_any_on     flag_diabetes_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 559380 |   6804 | 566184
                 |  97.94 |   1.19 |  99.13
                 |  98.80 |   1.20 |
                 |  99.14 |  98.39 |
        ---------+--------+--------+
               1 |   4851 |    111 |   4962
                 |   0.85 |   0.02 |   0.87
                 |  97.76 |   2.24 |
                 |   0.86 |   1.61 |
        ---------+--------+--------+
        Total      564231     6915   571146
                    98.79     1.21   100.00

                   The SAS System         11:48 Mo


Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     44.0779    <.0001
  Likelihood Ratio Chi-Square    1     35.3568    <.0001
  Continuity Adj. Chi-Square     1     43.2166    <.0001
  Mantel-Haenszel Chi-Square     1     44.0778    <.0001
  Phi Coefficient                       0.0088
  Contingency Coefficient               0.0088
  Cramer's V                            0.0088


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    559380
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 571146

*/

/* Proportion of autoimmune DE ES events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_all_on=1 and flag_immuno_gene=1;
   tables flag_event_de*flag_diabetes_gene / chisq;
run;

/*
 
   Table of flag_event_de by flag_diabetes_gene

       flag_event_de
                 flag_diabetes_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |   1531 |     22 |   1553
                |  41.50 |   0.60 |  42.10
                |  98.58 |   1.42 |
                |  42.25 |  33.85 |
       ---------+--------+--------+
              1 |   2093 |     43 |   2136
                |  56.74 |   1.17 |  57.90
                |  97.99 |   2.01 |
                |  57.75 |  66.15 |
       ---------+--------+--------+
       Total        3624       65     3689
                   98.24     1.76   100.00

                  The SAS System         11:48 Monday

Statistics for Table of flag_event_de by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      1.8484    0.1740
   Likelihood Ratio Chi-Square    1      1.8908    0.1691
   Continuity Adj. Chi-Square     1      1.5198    0.2176
   Mantel-Haenszel Chi-Square     1      1.8479    0.1740
   Phi Coefficient                       0.0224
   Contingency Coefficient               0.0224
   Cramer's V                            0.0224


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      1531
             Left-sided Pr <= F          0.9328
             Right-sided Pr >= F         0.1080

             Table Probability (P)       0.0408
             Two-sided Pr <= P           0.2051

                     Sample Size = 3689

*/


/* Proportion of autoimmune DD_1 ES events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_1*flag_diabetes_gene / chisq;
run;

/*
 Table of flag_event_dd_1 by flag_diabetes_gene

      flag_event_dd_1
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   1531 |     22 |   1553
               |  35.63 |   0.51 |  36.14
               |  98.58 |   1.42 |
               |  36.40 |  24.18 |
      ---------+--------+--------+
             1 |   2675 |     69 |   2744
               |  62.25 |   1.61 |  63.86
               |  97.49 |   2.51 |
               |  63.60 |  75.82 |
      ---------+--------+--------+
      Total        4206       91     4297
                  97.88     2.12   100.00

            Frequency Missing = 665

  Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1      5.7675    0.0163
      Likelihood Ratio Chi-Square    1      6.1313    0.0133
      Continuity Adj. Chi-Square     1      5.2500    0.0219
      Mantel-Haenszel Chi-Square     1      5.7661    0.0163
      Phi Coefficient                       0.0366
      Contingency Coefficient               0.0366
      Cramer's V                            0.0366


 Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      1531
               Left-sided Pr <= F          0.9951
               Right-sided Pr >= F         0.0095

               Table Probability (P)       0.0046
               Two-sided Pr <= P           0.0153

                  Effective Sample Size = 4297

*/


/* Proportion of autoimmune DD_2 ES events that are diabetes events */


proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq;
run;


/*
 Table of flag_event_dd_2 by flag_diabetes_gene

      flag_event_dd_2
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   1531 |     22 |   1553
               |  30.85 |   0.44 |  31.30
               |  98.58 |   1.42 |
               |  31.56 |  19.82 |
      ---------+--------+--------+
             1 |   3320 |     89 |   3409
               |  66.91 |   1.79 |  68.70
               |  97.39 |   2.61 |
               |  68.44 |  80.18 |
      ---------+--------+--------+
      Total        4851      111     4962
                  97.76     2.24   100.00


  Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1      6.9566    0.0084
      Likelihood Ratio Chi-Square    1      7.5535    0.0060
      Continuity Adj. Chi-Square     1      6.4213    0.0113
      Mantel-Haenszel Chi-Square     1      6.9552    0.0084
      Phi Coefficient                       0.0374
      Contingency Coefficient               0.0374
      Cramer's V                            0.0374


                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)      1531
                Left-sided Pr <= F          0.9978
                Right-sided Pr >= F         0.0044

                Table Probability (P)       0.0022
                Two-sided Pr <= P           0.0093

                        Sample Size = 4962


*/





/* Proportion of detected autoimmune IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_any_on*flag_diabetes_gene / chisq;
run;


/*
Table of flag_any_on by flag_diabetes_gene

   flag_any_on     flag_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  15919 |    261 |  16180
            |  76.26 |   1.25 |  77.51
            |  98.39 |   1.61 |
            |  77.65 |  69.97 |
   ---------+--------+--------+
          1 |   4582 |    112 |   4694
            |  21.95 |   0.54 |  22.49
            |  97.61 |   2.39 |
            |  22.35 |  30.03 |
   ---------+--------+--------+
   Total       20501      373    20874
               98.21     1.79   100.00


Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     12.3855    0.0004
  Likelihood Ratio Chi-Square    1     11.5729    0.0007
  Continuity Adj. Chi-Square     1     11.9490    0.0005
  Mantel-Haenszel Chi-Square     1     12.3849    0.0004
  Phi Coefficient                       0.0244
  Contingency Coefficient               0.0244
  Cramer's V                            0.0244


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     15919
            Left-sided Pr <= F          0.9997
            Right-sided Pr >= F         0.0004

            Table Probability (P)       0.0001
            Two-sided Pr <= P           0.0007

                   Sample Size = 20874


*/

/* Proportion of autoimmune DE IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_all_on=1 and flag_immuno_gene=1;
   tables flag_event_de*flag_diabetes_gene / chisq;
run;

/*
 Table of flag_event_de by flag_diabetes_gene

     flag_event_de
               flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   1883 |     23 |   1906
              |  61.44 |   0.75 |  62.19
              |  98.79 |   1.21 |
              |  62.39 |  48.94 |
     ---------+--------+--------+
            1 |   1135 |     24 |   1159
              |  37.03 |   0.78 |  37.81
              |  97.93 |   2.07 |
              |  37.61 |  51.06 |
     ---------+--------+--------+
     Total        3018       47     3065
                 98.47     1.53   100.00

 Statistics for Table of flag_event_de by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      3.5635    0.0591
    Likelihood Ratio Chi-Square    1      3.4515    0.0632
    Continuity Adj. Chi-Square     1      3.0143    0.0825
    Mantel-Haenszel Chi-Square     1      3.5624    0.0591
    Phi Coefficient                       0.0341
    Contingency Coefficient               0.0341
    Cramer's V                            0.0341


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      1883
              Left-sided Pr <= F          0.9780
              Right-sided Pr >= F         0.0428

              Table Probability (P)       0.0208
              Two-sided Pr <= P           0.0687

                      Sample Size = 3065


*/



/* Proportion of autoimmune DD_1 IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_1*flag_diabetes_gene / chisq;
run;

/*

 Table of flag_event_dd_1 by flag_diabetes_ge

      flag_event_dd_1
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   1883 |     23 |   1906
               |  48.15 |   0.59 |  48.73
               |  98.79 |   1.21 |
               |  49.15 |  28.75 |
      ---------+--------+--------+
             1 |   1948 |     57 |   2005
               |  49.81 |   1.46 |  51.27
               |  97.16 |   2.84 |
               |  50.85 |  71.25 |
      ---------+--------+--------+
      Total        3831       80     3911
                  97.95     2.05   100.00

 Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     13.0552    0.0003
     Likelihood Ratio Chi-Square    1     13.5165    0.0002
     Continuity Adj. Chi-Square     1     12.2514    0.0005
     Mantel-Haenszel Chi-Square     1     13.0519    0.0003
     Phi Coefficient                       0.0578
     Contingency Coefficient               0.0577
     Cramer's V                            0.0578

   Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)      1883
                 Left-sided Pr <= F          0.9999
                 Right-sided Pr >= F         0.0002

                 Table Probability (P)       0.0001
                 Two-sided Pr <= P           0.0003

                    Effective Sample Size = 3911

*/


/* Proportion of autoimmune DD_2 IR events that are diabetes events */


proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq;
run;


/*
 Table of flag_event_dd_2 by flag_diabetes_gene

      flag_event_dd_2
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   1883 |     23 |   1906
               |  40.12 |   0.49 |  40.61
               |  98.79 |   1.21 |
               |  41.10 |  20.54 |
      ---------+--------+--------+
             1 |   2699 |     89 |   2788
               |  57.50 |   1.90 |  59.39
               |  96.81 |   3.19 |
               |  58.90 |  79.46 |
      ---------+--------+--------+
      Total        4582      112     4694
                  97.61     2.39   100.00

 Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     19.1620    <.0001
     Likelihood Ratio Chi-Square    1     20.9127    <.0001
     Continuity Adj. Chi-Square     1     18.3190    <.0001
     Mantel-Haenszel Chi-Square     1     19.1579    <.0001
     Phi Coefficient                       0.0639
     Contingency Coefficient               0.0638
     Cramer's V                            0.0639


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      1883
               Left-sided Pr <= F          1.0000
               Right-sided Pr >= F         <.0001

               Table Probability (P)       <.0001
               Two-sided Pr <= P           <.0001

                       Sample Size = 4694

 

*/

/* Among detected junctions, what is the proportion of unannotated from autoimmune genes vs all genes? */


proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=0;
   tables flag_junction_annotated*flag_immuno_gene / chisq ;
run;

/*
  Table of flag_junction_annotated by flag_immuno_gene

         flag_junction_annotated
                   flag_immuno_gene

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |  34262 |   4455 |  38717
                  |  20.27 |   2.64 |  22.91
                  |  88.49 |  11.51 |
                  |  22.99 |  22.33 |
         ---------+--------+--------+
                1 | 114783 |  15495 | 130278
                  |  67.92 |   9.17 |  77.09
                  |  88.11 |  11.89 |
                  |  77.01 |  77.67 |
         ---------+--------+--------+
         Total      149045    19950   168995
                     88.19    11.81   100.00

 Statistics for Table of flag_junction_annotated by flag_immuno_gene

        Statistic                     DF       Value      Prob
        ------------------------------------------------------
        Chi-Square                     1      4.2984    0.0381
        Likelihood Ratio Chi-Square    1      4.3208    0.0376
        Continuity Adj. Chi-Square     1      4.2613    0.0390
        Mantel-Haenszel Chi-Square     1      4.2984    0.0381
        Phi Coefficient                       0.0050
        Contingency Coefficient               0.0050
        Cramer's V                            0.0050


                         Fisher's Exact Test
                  ----------------------------------
                  Cell (1,1) Frequency (F)     34262
                  Left-sided Pr <= F          0.9815
                  Right-sided Pr >= F         0.0193

                  Table Probability (P)       0.0008
                  Two-sided Pr <= P           0.0383

                         Sample Size = 168995


*/

/* Among detected junctions in autoimmune genes, what is the proportion of unannotated from diabetes genes ?*/

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_junction_annotated*flag_diabetes_gene / chisq ;
run;

/*

Table of flag_junction_annotated by flag_diabetes_gene

         flag_junction_annotated
                   flag_diabetes_gene

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |   4369 |     86 |   4455
                  |  21.90 |   0.43 |  22.33
                  |  98.07 |   1.93 |
                  |  22.29 |  24.43 |
         ---------+--------+--------+
                1 |  15229 |    266 |  15495
                  |  76.34 |   1.33 |  77.67
                  |  98.28 |   1.72 |
                  |  77.71 |  75.57 |
         ---------+--------+--------+
         Total       19598      352    19950
                     98.24     1.76   100.00


 Statistics for Table of flag_junction_annotated by flag_diabetes_g

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      0.9119    0.3396
         Likelihood Ratio Chi-Square    1      0.8929    0.3447
         Continuity Adj. Chi-Square     1      0.7928    0.3733
         Mantel-Haenszel Chi-Square     1      0.9119    0.3396
         Phi Coefficient                      -0.0068
         Contingency Coefficient               0.0068
         Cramer's V                           -0.0068


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)      4369
                   Left-sided Pr <= F          0.1860
                   Right-sided Pr >= F         0.8460

                   Table Probability (P)       0.0319
                   Two-sided Pr <= P           0.3333

                          Sample Size = 19950



*/


/* Among detected events, what is the proportion of IR events from autoimmune genes vs all genes? */


proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_intron_retention*flag_immuno_gene / chisq ;
run;

/*
Table of flag_intron_retention by flag_immuno_gene

       flag_intron_retention
                 flag_immuno_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 | 149045 |  19950 | 168995
                |  71.01 |   9.51 |  80.52
                |  88.19 |  11.81 |
                |  80.46 |  80.95 |
       ---------+--------+--------+
              1 |  36198 |   4694 |  40892
                |  17.25 |   2.24 |  19.48
                |  88.52 |  11.48 |
                |  19.54 |  19.05 |
       ---------+--------+--------+
       Total      185243    24644   209887
                   88.26    11.74   100.00

 Statistics for Table of flag_intron_retention by flag_immuno_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      3.3780    0.0661
       Likelihood Ratio Chi-Square    1      3.3947    0.0654
       Continuity Adj. Chi-Square     1      3.3466    0.0673
       Mantel-Haenszel Chi-Square     1      3.3779    0.0661
       Phi Coefficient                      -0.0040
       Contingency Coefficient               0.0040
       Cramer's V                           -0.0040


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)    149045
                 Left-sided Pr <= F          0.0334
                 Right-sided Pr >= F         0.9678

                 Table Probability (P)       0.0013
                 Two-sided Pr <= P           0.0670

                        Sample Size = 209887

*/

/* Among detected events in autoimmune genes, what is the proportion of IR events from diabetes genes ?*/


proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_intron_retention*flag_diabetes_gene / chisq ;
run;

/*
  Table of flag_intron_retention by flag_diabetes_gene

          flag_intron_retention
                    flag_diabetes_gene

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 |  19598 |    352 |  19950
                   |  79.52 |   1.43 |  80.95
                   |  98.24 |   1.76 |
                   |  81.05 |  75.86 |
          ---------+--------+--------+
                 1 |   4582 |    112 |   4694
                   |  18.59 |   0.45 |  19.05
                   |  97.61 |   2.39 |
                   |  18.95 |  24.14 |
          ---------+--------+--------+
          Total       24180      464    24644
                      98.12     1.88   100.00

 Statistics for Table of flag_intron_retention by flag_diabetes_gene

        Statistic                     DF       Value      Prob
        ------------------------------------------------------
        Chi-Square                     1      7.9481    0.0048
        Likelihood Ratio Chi-Square    1      7.4829    0.0062
        Continuity Adj. Chi-Square     1      7.6152    0.0058
        Mantel-Haenszel Chi-Square     1      7.9478    0.0048
        Phi Coefficient                       0.0180
        Contingency Coefficient               0.0180
        Cramer's V                            0.0180


                         Fisher's Exact Test
                  ----------------------------------
                  Cell (1,1) Frequency (F)     19598
                  Left-sided Pr <= F          0.9975
                  Right-sided Pr >= F         0.0036

                  Table Probability (P)       0.0010
                  Two-sided Pr <= P           0.0060

                         Sample Size = 24644




*/


/* Proportion of DE unannotated junctions among all genes that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_all_on=1;
   tables flag_event_de*flag_immuno_gene / chisq;
run;

/*

 Table of flag_event_de by flag_immuno_gene

    flag_event_de     flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  12359 |   1400 |  13759
             |  42.11 |   4.77 |  46.88
             |  89.82 |  10.18 |
             |  47.56 |  41.62 |
    ---------+--------+--------+
           1 |  13627 |   1964 |  15591
             |  46.43 |   6.69 |  53.12
             |  87.40 |  12.60 |
             |  52.44 |  58.38 |
    ---------+--------+--------+
    Total       25986     3364    29350
                88.54    11.46   100.00

               The SAS System         11:48

 Statistics for Table of flag_event_de by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     42.2444    <.0001
   Likelihood Ratio Chi-Square    1     42.4781    <.0001
   Continuity Adj. Chi-Square     1     42.0061    <.0001
   Mantel-Haenszel Chi-Square     1     42.2430    <.0001
   Phi Coefficient                       0.0379
   Contingency Coefficient               0.0379
   Cramer's V                            0.0379


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     12359
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 29350



*/

/* Proportion of DD_1 unannotated junctions among all genes that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_any_on=1;
   tables flag_event_dd_1*flag_immuno_gene / chisq;
run;


/*

  Table of flag_event_dd_1 by flag_immuno_gene

     flag_event_dd_1
               flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  12359 |   1400 |  13759
              |  36.37 |   4.12 |  40.49
              |  89.82 |  10.18 |
              |  41.06 |  36.06 |
     ---------+--------+--------+
            1 |  17743 |   2482 |  20225
              |  52.21 |   7.30 |  59.51
              |  87.73 |  12.27 |
              |  58.94 |  63.94 |
     ---------+--------+--------+
     Total       30102     3882    33984
                 88.58    11.42   100.00.

 Statistics for Table of flag_event_dd_1 by flag_immuno_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     35.5800    <.0001
    Likelihood Ratio Chi-Square    1     35.9921    <.0001
    Continuity Adj. Chi-Square     1     35.3731    <.0001
    Mantel-Haenszel Chi-Square     1     35.5789    <.0001
    Phi Coefficient                       0.0324
    Contingency Coefficient               0.0323
    Cramer's V                            0.0324

 Statistics for Table of flag_event_dd_1 by flag_immuno_gene

                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     12359
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                Effective Sample Size = 33984

*/

/* Proportion of DD_2 unannotated junctions among all genes that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;

/*
 Table of flag_event_dd_2 by flag_immuno_gene

     flag_event_dd_2
               flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  12359 |   1400 |  13759
              |  31.92 |   3.62 |  35.54
              |  89.82 |  10.18 |
              |  36.07 |  31.43 |
     ---------+--------+--------+
            1 |  21903 |   3055 |  24958
              |  56.57 |   7.89 |  64.46
              |  87.76 |  12.24 |
              |  63.93 |  68.57 |
     ---------+--------+--------+
     Total       34262     4455    38717
                 88.49    11.51   100.00

 
 Statistics for Table of flag_event_dd_2 by flag_immuno_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     37.1577    <.0001
    Likelihood Ratio Chi-Square    1     37.7758    <.0001
    Continuity Adj. Chi-Square     1     36.9551    <.0001
    Mantel-Haenszel Chi-Square     1     37.1567    <.0001
    Phi Coefficient                       0.0310
    Contingency Coefficient               0.0310
    Cramer's V                            0.0310


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     12359
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 38717

*/

/* Proportion of DE unannotated junctions among autoimmune genes that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_all_on=1 and flag_immuno_gene=1;
   tables flag_event_de*flag_diabetes_gene / chisq;
run;


/*
 Table of flag_event_de by flag_diabetes_gene

     flag_event_de
               flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   1387 |     13 |   1400
              |  41.23 |   0.39 |  41.62
              |  99.07 |   0.93 |
              |  41.75 |  30.95 |
     ---------+--------+--------+
            1 |   1935 |     29 |   1964
              |  57.52 |   0.86 |  58.38
              |  98.52 |   1.48 |
              |  58.25 |  69.05 |
     ---------+--------+--------+
     Total        3322       42     3364
                 98.75     1.25   100.00


 Statistics for Table of flag_event_de by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      1.9909    0.1582
    Likelihood Ratio Chi-Square    1      2.0581    0.1514
    Continuity Adj. Chi-Square     1      1.5712    0.2100
    Mantel-Haenszel Chi-Square     1      1.9903    0.1583
    Phi Coefficient                       0.0243
    Contingency Coefficient               0.0243
    Cramer's V                            0.0243


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      1387
              Left-sided Pr <= F          0.9437
              Right-sided Pr >= F         0.1038

              Table Probability (P)       0.0475
              Two-sided Pr <= P           0.2071

                      Sample Size = 3364

*/

/* Proportion of DD_1 unannotated junctions among autoimmune genes that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_1*flag_diabetes_gene / chisq;
run;

/*

     Table of flag_event_dd_1 by flag_diabetes_gene

          flag_event_dd_1
                    flag_diabetes_gene

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 |   1387 |     13 |   1400
                   |  35.73 |   0.33 |  36.06
                   |  99.07 |   0.93 |
                   |  36.33 |  20.31 |
          ---------+--------+--------+
                 1 |   2431 |     51 |   2482
                   |  62.62 |   1.31 |  63.94
                   |  97.95 |   2.05 |
                   |  63.67 |  79.69 |
          ---------+--------+--------+
          Total        3818       64     3882
                      98.35     1.65   100.00

  Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1      7.0019    0.0081
      Likelihood Ratio Chi-Square    1      7.6532    0.0057
      Continuity Adj. Chi-Square     1      6.3246    0.0119
      Mantel-Haenszel Chi-Square     1      7.0001    0.0082
      Phi Coefficient                       0.0425
      Contingency Coefficient               0.0424
      Cramer's V                            0.0425

  Statistics for Table of flag_event_dd_1 by flag_diabetes_gene

                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)      1387
                Left-sided Pr <= F          0.9981
                Right-sided Pr >= F         0.0046

                Table Probability (P)       0.0027
                Two-sided Pr <= P           0.0082

                   Effective Sample Size = 3882


*/

/* Proportion of DD_2 unannotated junctions among autoimmune genes that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq;
run;

/*
   Table of flag_event_dd_2 by flag_diabetes_gene

      flag_event_dd_2
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   1387 |     13 |   1400
               |  31.13 |   0.29 |  31.43
               |  99.07 |   0.93 |
               |  31.75 |  15.12 |
      ---------+--------+--------+
             1 |   2982 |     73 |   3055
               |  66.94 |   1.64 |  68.57
               |  97.61 |   2.39 |
               |  68.25 |  84.88 |
      ---------+--------+--------+
      Total        4369       86     4455
                  98.07     1.93   100.00

 Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     10.8238    0.0010
     Likelihood Ratio Chi-Square    1     12.3311    0.0004
     Continuity Adj. Chi-Square     1     10.0658    0.0015
     Mantel-Haenszel Chi-Square     1     10.8214    0.0010
     Phi Coefficient                       0.0493
     Contingency Coefficient               0.0492
     Cramer's V                            0.0493


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      1387
               Left-sided Pr <= F          0.9999
               Right-sided Pr >= F         0.0004

               Table Probability (P)       0.0003
               Two-sided Pr <= P           0.0006


*/



/******* DETECTED AMONG ALL ********/

/* Proportion of exon skipping of all detected events */

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_exonskip*flag_immuno_gene / chisq;
run;

/*
   Table of flag_exonskip by flag_immuno_gene

     flag_exonskip     flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 | 144124 |  19682 | 163806
              |  68.67 |   9.38 |  78.04
              |  87.98 |  12.02 |
              |  77.80 |  79.87 |
     ---------+--------+--------+
            1 |  41119 |   4962 |  46081
              |  19.59 |   2.36 |  21.96
              |  89.23 |  10.77 |
              |  22.20 |  20.13 |
     ---------+--------+--------+
     Total      185243    24644   209887
                 88.26    11.74   100.00


 Statistics for Table of flag_exonskip by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     54.0035    <.0001
   Likelihood Ratio Chi-Square    1     54.9794    <.0001
   Continuity Adj. Chi-Square     1     53.8832    <.0001
   Mantel-Haenszel Chi-Square     1     54.0032    <.0001
   Phi Coefficient                      -0.0160
   Contingency Coefficient               0.0160
   Cramer's V                           -0.0160


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    144124
             Left-sided Pr <= F          <.0001
             Right-sided Pr >= F         1.0000

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 209887


*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_exonskip*flag_diabetes_gene / chisq;
run;


/*
   Table of flag_exonskip by flag_diabetes_gene

      flag_exonskip
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  19329 |    353 |  19682
               |  78.43 |   1.43 |  79.87
               |  98.21 |   1.79 |
               |  79.94 |  76.08 |
      ---------+--------+--------+
             1 |   4851 |    111 |   4962
               |  19.68 |   0.45 |  20.13
               |  97.76 |   2.24 |
               |  20.06 |  23.92 |
      ---------+--------+--------+
      Total       24180      464    24644
                  98.12     1.88   100.00

                 The SAS System         11:48 Mo

Statistics for Table of flag_exonskip by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      4.2191    0.0400
   Likelihood Ratio Chi-Square    1      4.0429    0.0444
   Continuity Adj. Chi-Square     1      3.9824    0.0460
   Mantel-Haenszel Chi-Square     1      4.2189    0.0400
   Phi Coefficient                       0.0131
   Contingency Coefficient               0.0131
   Cramer's V                            0.0131


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     19329
             Left-sided Pr <= F          0.9811
             Right-sided Pr >= F         0.0247

             Table Probability (P)       0.0058
             Two-sided Pr <= P           0.0467

                    Sample Size = 24644

*/



/* Proportion of exon skipping of all detected junctions */

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=0;
   tables flag_exonskip*flag_immuno_gene / chisq;
run;

/*
  Table of flag_exonskip by flag_immuno_gene

     flag_exonskip     flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 | 107926 |  14988 | 122914
              |  63.86 |   8.87 |  72.73
              |  87.81 |  12.19 |
              |  72.41 |  75.13 |
     ---------+--------+--------+
            1 |  41119 |   4962 |  46081
              |  24.33 |   2.94 |  27.27
              |  89.23 |  10.77 |
              |  27.59 |  24.87 |
     ---------+--------+--------+
     Total      149045    19950   168995
                 88.19    11.81   100.00

    Statistics for Table of flag_exonskip by flag_immuno_gene

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1     65.4504    <.0001
      Likelihood Ratio Chi-Square    1     66.5481    <.0001
      Continuity Adj. Chi-Square     1     65.3136    <.0001
      Mantel-Haenszel Chi-Square     1     65.4501    <.0001
      Phi Coefficient                      -0.0197
      Contingency Coefficient               0.0197
      Cramer's V                           -0.0197


                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)    107926
                Left-sided Pr <= F          <.0001
                Right-sided Pr >= F         1.0000

                Table Probability (P)       <.0001
                Two-sided Pr <= P           <.0001

                       Sample Size = 168995

*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_exonskip*flag_diabetes_gene / chisq;
run;


/*
 Table of flag_exonskip by flag_diabetes_gene

     flag_exonskip
               flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  14747 |    241 |  14988
              |  73.92 |   1.21 |  75.13
              |  98.39 |   1.61 |
              |  75.25 |  68.47 |
     ---------+--------+--------+
            1 |   4851 |    111 |   4962
              |  24.32 |   0.56 |  24.87
              |  97.76 |   2.24 |
              |  24.75 |  31.53 |
     ---------+--------+--------+
     Total       19598      352    19950
                 98.24     1.76   100.00

                The SAS System         11:48 Mo

 Statistics for Table of flag_exonskip by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      8.5105    0.0035
    Likelihood Ratio Chi-Square    1      8.0791    0.0045
    Continuity Adj. Chi-Square     1      8.1515    0.0043
    Mantel-Haenszel Chi-Square     1      8.5101    0.0035
    Phi Coefficient                       0.0207
    Contingency Coefficient               0.0206
    Cramer's V                            0.0207


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     14747
              Left-sided Pr <= F          0.9982
              Right-sided Pr >= F         0.0026

              Table Probability (P)       0.0008
              Two-sided Pr <= P           0.0042

                     Sample Size = 19950

*/

/* Proportion of IR events of all detected events */

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_intron_retention*flag_immuno_gene / chisq;
run;


/*
  Table of flag_intron_retention by flag_immuno_gene

        flag_intron_retention
                  flag_immuno_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 149045 |  19950 | 168995
                 |  71.01 |   9.51 |  80.52
                 |  88.19 |  11.81 |
                 |  80.46 |  80.95 |
        ---------+--------+--------+
               1 |  36198 |   4694 |  40892
                 |  17.25 |   2.24 |  19.48
                 |  88.52 |  11.48 |
                 |  19.54 |  19.05 |
        ---------+--------+--------+
        Total      185243    24644   209887
                    88.26    11.74   100.00

                   The SAS System         11:48 Mond

                 The FREQ Procedure

 Statistics for Table of flag_intron_retention by flag_immuno_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      3.3780    0.0661
       Likelihood Ratio Chi-Square    1      3.3947    0.0654
       Continuity Adj. Chi-Square     1      3.3466    0.0673
       Mantel-Haenszel Chi-Square     1      3.3779    0.0661
       Phi Coefficient                      -0.0040
       Contingency Coefficient               0.0040
       Cramer's V                           -0.0040


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)    149045
                 Left-sided Pr <= F          0.0334
                 Right-sided Pr >= F         0.9678

                 Table Probability (P)       0.0013
                 Two-sided Pr <= P           0.0670

                        Sample Size = 209887


*/


* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_intron_retention*flag_diabetes_gene / chisq;
run;


/*
 
Table of flag_intron_retention by flag_diabetes_gene

        flag_intron_retention
                  flag_diabetes_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |  19598 |    352 |  19950
                 |  79.52 |   1.43 |  80.95
                 |  98.24 |   1.76 |
                 |  81.05 |  75.86 |
        ---------+--------+--------+
               1 |   4582 |    112 |   4694
                 |  18.59 |   0.45 |  19.05
                 |  97.61 |   2.39 |
                 |  18.95 |  24.14 |
        ---------+--------+--------+
        Total       24180      464    24644
                    98.12     1.88   100.00

                   The SAS System         11:48 Monday,

Statistics for Table of flag_intron_retention by flag_diabetes_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      7.9481    0.0048
       Likelihood Ratio Chi-Square    1      7.4829    0.0062
       Continuity Adj. Chi-Square     1      7.6152    0.0058
       Mantel-Haenszel Chi-Square     1      7.9478    0.0048
       Phi Coefficient                       0.0180
       Contingency Coefficient               0.0180
       Cramer's V                            0.0180


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)     19598
                 Left-sided Pr <= F          0.9975
                 Right-sided Pr >= F         0.0036

                 Table Probability (P)       0.0010
                 Two-sided Pr <= P           0.0060

                        Sample Size = 24644



*/

/* Proportion of Unannotated events of all detected events */


data splicing_clean_w_de_flags2;
   set splicing_clean_w_de_flags;
   if flag_junction_annotated=0 and flag_intron_retention=0 then flag_junction_unannotated=1;
   else flag_junction_unannotated=0;
run;

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1;
   tables flag_junction_unannotated*flag_immuno_gene / chisq ;
run;


/*

 Table of flag_junction_unannotated by flag_immuno_gene

          flag_junction_unannotated
                    flag_immuno_gene

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 | 150981 |  20189 | 171170
                   |  71.93 |   9.62 |  81.55
                   |  88.21 |  11.79 |
                   |  81.50 |  81.92 |
          ---------+--------+--------+
                 1 |  34262 |   4455 |  38717
                   |  16.32 |   2.12 |  18.45
                   |  88.49 |  11.51 |
                   |  18.50 |  18.08 |
          ---------+--------+--------+
          Total      185243    24644   209887
                      88.26    11.74   100.00


 Statistics for Table of flag_junction_unannotated by flag_immuno_gene

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      2.5296    0.1117
         Likelihood Ratio Chi-Square    1      2.5410    0.1109
         Continuity Adj. Chi-Square     1      2.5019    0.1137
         Mantel-Haenszel Chi-Square     1      2.5296    0.1117
         Phi Coefficient                      -0.0035
         Contingency Coefficient               0.0035
         Cramer's V                           -0.0035


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)    150981
                   Left-sided Pr <= F          0.0566
                   Right-sided Pr >= F         0.9454

                   Table Probability (P)       0.0020
                   Two-sided Pr <= P           0.1136

                          Sample Size = 209887

*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_junction_unannotated*flag_diabetes_gene / chisq ;
run;


/*
  Table of flag_junction_unannotated by flag_diabetes_gene

            flag_junction_unannotated
                      flag_diabetes_gene

            Frequency|
            Percent  |
            Row Pct  |
            Col Pct  |       0|       1|  Total
            ---------+--------+--------+
                   0 |  19811 |    378 |  20189
                     |  80.39 |   1.53 |  81.92
                     |  98.13 |   1.87 |
                     |  81.93 |  81.47 |
            ---------+--------+--------+
                   1 |   4369 |     86 |   4455
                     |  17.73 |   0.35 |  18.08
                     |  98.07 |   1.93 |
                     |  18.07 |  18.53 |
            ---------+--------+--------+
            Total       24180      464    24644
                        98.12     1.88   100.00

                       The SAS System         11:48 Monday
Statistics for Table of flag_junction_unannotated by flag_diabetes_gene

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      0.0667    0.7962
         Likelihood Ratio Chi-Square    1      0.0663    0.7968
         Continuity Adj. Chi-Square     1      0.0390    0.8435
         Mantel-Haenszel Chi-Square     1      0.0667    0.7962
         Phi Coefficient                       0.0016
         Contingency Coefficient               0.0016
         Cramer's V                            0.0016


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)     19811
                   Left-sided Pr <= F          0.6294
                   Right-sided Pr >= F         0.4171

                   Table Probability (P)       0.0465
                   Two-sided Pr <= P           0.8075

                          Sample Size = 24644


*/


/* IR+Unannot vs Annotated junctions */


* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1;
   tables flag_junction_annotated*flag_immuno_gene / chisq ;
run;

/*

 Table of flag_junction_annotated by flag_immuno_gene

         flag_junction_annotated
                   flag_immuno_gene

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |  70460 |   9149 |  79609
                  |  33.57 |   4.36 |  37.93
                  |  88.51 |  11.49 |
                  |  38.04 |  37.12 |
         ---------+--------+--------+
                1 | 114783 |  15495 | 130278
                  |  54.69 |   7.38 |  62.07
                  |  88.11 |  11.89 |
                  |  61.96 |  62.88 |
         ---------+--------+--------+
         Total      185243    24644   209887
                     88.26    11.74   100.00

Statistics for Table of flag_junction_annotated by flag_immuno_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      7.6819    0.0056
       Likelihood Ratio Chi-Square    1      7.7007    0.0055
       Continuity Adj. Chi-Square     1      7.6433    0.0057
       Mantel-Haenszel Chi-Square     1      7.6819    0.0056
       Phi Coefficient                       0.0060
       Contingency Coefficient               0.0060
       Cramer's V                            0.0060


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)     70460
                 Left-sided Pr <= F          0.9973
                 Right-sided Pr >= F         0.0028

                 Table Probability (P)       0.0001
                 Two-sided Pr <= P           0.0057

                        Sample Size = 209887

*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_junction_annotated*flag_diabetes_gene / chisq ;
run;


/*
  Table of flag_junction_annotated by flag_diabetes_gene

           flag_junction_annotated
                     flag_diabetes_gene

           Frequency|
           Percent  |
           Row Pct  |
           Col Pct  |       0|       1|  Total
           ---------+--------+--------+
                  0 |   8951 |    198 |   9149
                    |  36.32 |   0.80 |  37.12
                    |  97.84 |   2.16 |
                    |  37.02 |  42.67 |
           ---------+--------+--------+
                  1 |  15229 |    266 |  15495
                    |  61.80 |   1.08 |  62.88
                    |  98.28 |   1.72 |
                    |  62.98 |  57.33 |
           ---------+--------+--------+
           Total       24180      464    24644
                       98.12     1.88   100.00

                      The SAS System         11:48 Monday,

 Statistics for Table of flag_junction_annotated by flag_diabetes_gene

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      6.2354    0.0125
         Likelihood Ratio Chi-Square    1      6.1269    0.0133
         Continuity Adj. Chi-Square     1      5.9955    0.0143
         Mantel-Haenszel Chi-Square     1      6.2352    0.0125
         Phi Coefficient                      -0.0159
         Contingency Coefficient               0.0159
         Cramer's V                           -0.0159


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)      8951
                   Left-sided Pr <= F          0.0075
                   Right-sided Pr >= F         0.9942

                   Table Probability (P)       0.0018
                   Two-sided Pr <= P           0.0133

                          Sample Size = 24644


*/

