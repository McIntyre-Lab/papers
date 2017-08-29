/* Counting splicing events for T1D genes */

/* Want counts for total, exon-skipping, alt donor, alt acceptor, alt donor+alt acceptor, intron retention, annotated junc, unannotated junc */
/* Want for each cell type, and overlaps */
/* And for events shared in all three, which are DE */

libname splicing '/mnt/data/splicing/';
libname con '/home/jrbnewman/concannon/sas_data';

data immunogenes;
   set con.immunogene_flags;
   if flag_diabetes_gene=1;
   keep gene_id;
run;

proc sort data=immunogenes;
   by gene_id;
proc sort data=splicing.splicing_results_w_annot_fdr;
   by gene_id;
run;


data splicing.splicing_results_annot_t1d;
    merge immunogenes (in=in1) splicing.splicing_results_w_annot_fdr (in=in2);
    by gene_id;
    if in1 and in2 then output;
run;

ods listing; ods html close;

/* Counts for all events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_all;
run;

/* Counts for ES events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_exonskip=1;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_es;
run;

/* Counts for AD-only events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_ad;
run;

/* Counts for AA-only events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_aa;
run;


/* Counts for AD+AA events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_ad_aa;
run;

/* Counts for IR events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_intron_retention=1;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_ir;
run;

/* Counts for annotated junc events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_junction_annotated=1;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_annot;
run;

/* Counts for unannotated junc (no IR) events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_junction_annotated=0 and flag_intron_retention=0;
    tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=cell_counts_unannot;
run;

/* Counts for DE events */
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1;
    tables flag_anova_fdr_05 / out=cell_counts_all_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_exonskip=1;
    tables flag_anova_fdr_05 / out=cell_counts_es_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0;
    tables flag_anova_fdr_05 / out=cell_counts_ad_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0;
    tables flag_anova_fdr_05 / out=cell_counts_aa_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0;
    tables flag_anova_fdr_05 / out=cell_counts_ad_aa_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_intron_retention=1;
    tables flag_anova_fdr_05 / out=cell_counts_ir_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_junction_annotated=1;
    tables flag_anova_fdr_05 / out=cell_counts_annot_sig;
proc freq data=splicing.splicing_results_annot_t1d noprint;
    where flag_all_on=1 and flag_junction_annotated=0 and flag_intron_retention=0;
    tables flag_anova_fdr_05 / out=cell_counts_unannot_sig;
run;


/* Stack all counts together */

data cell_counts_all2;
    length event_type $10.;
    set cell_counts_all;
    event_type='all';    
    drop percent;
run;

data cell_counts_es2;
    length event_type $10.;
    set cell_counts_es;
    event_type='ES';    
    drop percent;
run;

data cell_counts_ad2;
    length event_type $10.;
    set cell_counts_ad;
    event_type='AD';    
    drop percent;
run;

data cell_counts_aa2;
    length event_type $10.;
    set cell_counts_aa;
    event_type='AA';    
    drop percent;
run;

data cell_counts_ad_aa2;
    length event_type $10.;
    set cell_counts_ad_aa;
    event_type='AD_AA';    
    drop percent;
run;

data cell_counts_ir2;
    length event_type $10.;
    set cell_counts_ir;
    event_type='IR';    
    drop percent;
run;

data cell_counts_annot2;
    length event_type $10.;
    set cell_counts_annot;
    event_type='Annot';    
    drop percent;
run;

data cell_counts_unannot2;
    length event_type $10.;
    set cell_counts_unannot;
    event_type='Unannot';    
    drop percent;
run;


data all_event_counts;
   set cell_counts_all2 cell_counts_es2 cell_counts_ad2 cell_counts_aa2
   cell_counts_ad_aa2 cell_counts_ir2 cell_counts_annot2 cell_counts_unannot2;
run;


/* Stack significant counts */
data cell_counts_all_sig2;
   set cell_counts_all_sig;
   length event_type $10.;
   event_type='all';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_es_sig2;
   set cell_counts_es_sig;
   length event_type $10.;
   event_type='ES';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_ad_sig2;
   set cell_counts_ad_sig;
   length event_type $10.;
   event_type='AD';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_aa_sig2;
   set cell_counts_aa_sig;
   length event_type $10.;
   event_type='AA';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_ad_aa_sig2;
   set cell_counts_ad_aa_sig;
   length event_type $10.;
   event_type='AD_AA';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_ir_sig2;
   set cell_counts_ir_sig;
   length event_type $10.;
   event_type='IR';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_annot_sig2;
   set cell_counts_annot_sig;
   length event_type $10.;
   event_type='Annot';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data cell_counts_unannot_sig2;
   set cell_counts_unannot_sig;
   length event_type $10.;
   event_type='Unannot';
   if flag_anova_fdr_05 ne 1 then delete;
   rename count=num_sig_events;
   drop flag_anova_fdr_05 percent;
run;

data sig_event_counts;
   set cell_counts_all_sig2 cell_counts_es_sig2 cell_counts_ad_sig2 cell_counts_aa_sig2
   cell_counts_ad_aa_sig2 cell_counts_ir_sig2 cell_counts_annot_sig2 cell_counts_unannot_sig2;
   flag_CD4_on=1;
   flag_CD8_on=1;
   flag_CD19_on=1;
run;


/* Merge counts */

proc sort data=all_event_counts;
   by flag_CD4_on flag_CD8_on flag_CD19_on event_type;
proc sort data=sig_event_counts;
   by flag_CD4_on flag_CD8_on flag_CD19_on event_type;
run;

data splicing.event_counts_summary_t1d oops;
   merge all_event_counts (in=in1) sig_event_counts;
   by flag_CD4_on flag_CD8_on flag_CD19_on event_type;
   if in1 then output splicing.event_counts_summary_t1d;
   else output oops;
run;


