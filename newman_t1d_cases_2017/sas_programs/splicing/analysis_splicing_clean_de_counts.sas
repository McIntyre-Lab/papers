/*** SPLICING DE COUNTS ***/

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';


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

