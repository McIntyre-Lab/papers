/* Get counts for fusions and genes */

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname splice '/mnt/data/splice';
libname splicing '/mnt/data/splicing';

/* Merge in Immunoflags */

data splicing_results;
   set splicing.splicing_results_w_annot_fdr;
   keep event_id flag_cd19_on flag_cd4_on flag_cd8_on flag_all_on
        flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05 gene_id
        flag_junction_annotated flag_intron_retention flag_exonskip flag_alt_donor flag_alt_acceptor
        flag_anova_fdr_05;
run;

data immunogene_flags;
   set con.immunogene_flags;
run;

proc sort data=immunogene_flags;
   by gene_id;
proc sort data=splicing_results;
   by gene_id;
run;

data splicing2gene_for_counts;
   merge immunogene_flags (in=in1) splicing_results (in=in2);
   by gene_id;
   if in1 and in2 then output;
   else if in1 then delete;
   else output;
run;


/* Counts for splicing events */

/* All genes */
proc freq data=splicing2gene_for_counts noprint;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_exonskip=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_es;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_ad;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_aa;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_alt_donor=1 and flag_alt_acceptor=1  and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_ada;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_junction_annotated=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_ea;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_eu;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_intron_retention=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_all_ir;
run;

/*
	ALL	ES	AD	AA	ADA	EA	EU	IR
CD4	4385	1244	431	404	554	1861	1082	2808
CD8	6514	1733	515	522	525	2544	1508	2462
CD19	10045	2461	894	943	1016	5187	2050	2808
CD4/8	16066	4269	1480	1476	1954	7485	3725	4856
CD4/19	1923	471	182	203	210	998	407	518
CD8/19	2991	723	245	262	319	1313	604	1074
ALL	170311	35203	17412	18076	24139	111051	29365	29895
OFF	6212508	5523465	1152798	1165223	813471	186015	5791501	234992
*/

/* Autoimmune genes */


proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_exonskip=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_es;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_ad;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_aa;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_alt_donor=1 and flag_alt_acceptor=1  and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_ada;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_junction_annotated=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_ea;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_eu;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_immuno_gene=1 and flag_intron_retention=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_ai_ir;
run;

/*
	ALL	ES	AD	AA	ADA	EA	EU	IR
CD4	483	157	43	51	71	176	137	170
CD8	676	141	50	52	56	270	115	291
CD19	1232	310	125	127	162	580	266	386
CD4/8	1936	539	197	173	304	894	457	585
CD4/19	219	53	17	24	26	95	56	68
CD8/19	334	73	27	29	43	143	60	131
ALL	19788	3689	2072	2179	3599	13337	3364	3087
OFF	638159	580142	128351	130417	110752	11927	608929	17303
*/


/* T1D genes */


proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_exonskip=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_es;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_ad;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_aa;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_alt_donor=1 and flag_alt_acceptor=1  and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_ada;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_junction_annotated=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_ea;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_eu;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_diabetes_gene=1 and flag_intron_retention=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=splicing_counts_on_t1d_ir;
run;

/*
	ALL	ES	AD	AA	ADA	EA	EU	IR
CD4	9	3	2	0	1	2	3	4
CD8	22	5	1	3	0	2	6	14
CD19	63	18	9	10	7	35	13	15
CD4/8	61	19	5	3	6	15	22	24
CD4/19	11	1	1	3	0	8	0	3
CD8/19	5	0	0	0	0	0	0	5
ALL	293	65	27	34	34	204	42	47
OFF	7713	6804	1536	173	1190	187	7265	261

*/


/* Counts for DE splicing */

/* All genes */
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1;
   tables flag_anova_fdr_05 / out=splicing_de_all;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_exonskip=1;
   tables flag_anova_fdr_05 / out=splicing_de_all_es;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
   tables flag_anova_fdr_05 / out=splicing_de_all_ad;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_anova_fdr_05 / out=splicing_de_all_aa;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0;
   tables flag_anova_fdr_05 / out=splicing_de_all_ada;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=1;
   tables flag_anova_fdr_05 / out=splicing_de_all_ea;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_anova_fdr_05 / out=splicing_de_all_eu;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=1;
   tables flag_anova_fdr_05 / out=splicing_de_all_ir;
run;


/* Autoimmune genes */

proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_es;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_ad;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_aa;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_ada;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=1 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_ea;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_eu;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_ai_ir;
run;

/* T1D genes */

proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_exonskip=1 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_es;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_ad;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_aa;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_ada;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=1 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_ea;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_eu;
run;
proc freq data=splicing2gene_for_counts noprint;
   where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=1 and flag_diabetes_gene=1;
   tables flag_anova_fdr_05 / out=splicing_de_t1d_ir;
run;


/*

	ALL	ES	AD	AA	ADA	EA	EU	IR
ALL	92266	18783	10080	10351	14188	65175	15598	11493
AI	11432	2136	1257	1365	2221	8300	1964	1168
T1D	205	43	20	28	20	152	29	24


*/


