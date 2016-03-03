/*** SPLICING COUNTS ***/

/* Redo counts and check proportions of DE and DD splicing events - cleaned dataset

Here I am going to look at the three DE/DD scenarios again:

1. DE events expressed in all three tissues
2. DE events expressed in all three tissues PLUS tissue-specific events
3. DE events expressed in all three tissues PLUS two-tissue events

* Need to make a makefile for all this!;
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


