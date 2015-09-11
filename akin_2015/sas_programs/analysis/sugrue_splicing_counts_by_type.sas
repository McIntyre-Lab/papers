/* Total events */

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=1;
   tables gene_id / out=sug_all_ej_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   tables gene_id / out=sug_all_eu_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=1;
   tables gene_id / out=sug_all_ir_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_exonskip=1;
   tables gene_id / out=sug_all_es_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=0;
   tables gene_id / out=sug_all_ad_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=0 and flag_alt_acceptor=1;
   tables gene_id / out=sug_all_aa_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=1;
   tables gene_id / out=sug_all_ada_by_gene;
run;

/* detected events */

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=1 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_ej_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_eu_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=1 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_ir_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_exonskip=1 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_es_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=0 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_ad_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=0 and flag_alt_acceptor=1 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_aa_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=1 and (flag_control_on=1 or flag_treat_on=1);
   tables gene_id / out=sug_dtct_ada_by_gene;
run;

/* Analyzed events */


proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=1 and flag_all_on=1;
   tables gene_id / out=sug_tested_ej_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0 and flag_all_on=1;
   tables gene_id / out=sug_tested_eu_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=1 and flag_all_on=1;
   tables gene_id / out=sug_tested_ir_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_exonskip=1 and flag_all_on=1;
   tables gene_id / out=sug_tested_es_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_all_on=1;
   tables gene_id / out=sug_tested_ad_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_all_on=1;
   tables gene_id / out=sug_tested_aa_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_all_on=1;
   tables gene_id / out=sug_tested_ada_by_gene;
run;

/* Significant events */



proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=1 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_ej_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_eu_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=1 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_ir_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_exonskip=1 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_es_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_ad_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_aa_by_gene;
run;

proc freq data=sgrloc.results_by_jnct_fdr noprint;
   where flag_intron_retention=0 and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_all_on=1 and flag_p05=1;
   tables gene_id / out=sug_de_ada_by_gene;
run;

data  sug_all_ej_by_gene2;
set sug_all_ej_by_gene;
   keep gene_id count;
   rename count=ej_total;
run;

data  sug_all_eu_by_gene2;
set sug_all_eu_by_gene;
   keep gene_id count;
   rename count=eu_total;
run;

data  sug_all_ir_by_gene2;
set sug_all_ir_by_gene;
   keep gene_id count;
   rename count=ir_total;
run;

data  sug_all_es_by_gene2;
set sug_all_es_by_gene;
   keep gene_id count;
   rename count=es_total;
run;

data  sug_all_ad_by_gene2;
set sug_all_ad_by_gene;
   keep gene_id count;
   rename count=ad_total;
run;

data  sug_all_aa_by_gene2;
set sug_all_aa_by_gene;
   keep gene_id count;
   rename count=aa_total;
run;

data  sug_all_ada_by_gene2;
set sug_all_ada_by_gene;
   keep gene_id count;
   rename count=ada_total;
run;



data  sug_dtct_ej_by_gene2;
set sug_dtct_ej_by_gene;
   keep gene_id count;
   rename count=ej_dtct;
run;

data  sug_dtct_eu_by_gene2;
set sug_dtct_eu_by_gene;
   keep gene_id count;
   rename count=eu_dtct;
run;

data  sug_dtct_ir_by_gene2;
set sug_dtct_ir_by_gene;
   keep gene_id count;
   rename count=ir_dtct;
run;

data  sug_dtct_es_by_gene2;
set sug_dtct_es_by_gene;
   keep gene_id count;
   rename count=es_dtct;
run;

data  sug_dtct_ad_by_gene2;
set sug_dtct_ad_by_gene;
   keep gene_id count;
   rename count=ad_dtct;
run;

data  sug_dtct_aa_by_gene2;
set sug_dtct_aa_by_gene;
   keep gene_id count;
   rename count=aa_dtct;
run;

data  sug_dtct_ada_by_gene2;
set sug_dtct_ada_by_gene;
   keep gene_id count;
   rename count=ada_dtct;
run;





data  sug_tested_ej_by_gene2;
set sug_tested_ej_by_gene;
   keep gene_id count;
   rename count=ej_tested;
run;

data  sug_tested_eu_by_gene2;
set sug_tested_eu_by_gene;
   keep gene_id count;
   rename count=eu_tested;
run;

data  sug_tested_ir_by_gene2;
set sug_tested_ir_by_gene;
   keep gene_id count;
   rename count=ir_tested;
run;

data  sug_tested_es_by_gene2;
set sug_tested_es_by_gene;
   keep gene_id count;
   rename count=es_tested;
run;

data  sug_tested_ad_by_gene2;
set sug_tested_ad_by_gene;
   keep gene_id count;
   rename count=ad_tested;
run;

data  sug_tested_aa_by_gene2;
set sug_tested_aa_by_gene;
   keep gene_id count;
   rename count=aa_tested;
run;

data  sug_tested_ada_by_gene2;
set sug_tested_ada_by_gene;
   keep gene_id count;
   rename count=ada_tested;
run;







data  sug_de_ej_by_gene2;
set sug_de_ej_by_gene;
   keep gene_id count;
   rename count=ej_de;
run;

data  sug_de_eu_by_gene2;
set sug_de_eu_by_gene;
   keep gene_id count;
   rename count=eu_de;
run;

data  sug_de_ir_by_gene2;
set sug_de_ir_by_gene;
   keep gene_id count;
   rename count=ir_de;
run;

data  sug_de_es_by_gene2;
set sug_de_es_by_gene;
   keep gene_id count;
   rename count=es_de;
run;

data  sug_de_ad_by_gene2;
set sug_de_ad_by_gene;
   keep gene_id count;
   rename count=ad_de;
run;

data  sug_de_aa_by_gene2;
set sug_de_aa_by_gene;
   keep gene_id count;
   rename count=aa_de;
run;

data  sug_de_ada_by_gene2;
set sug_de_ada_by_gene;
   keep gene_id count;
   rename count=ada_de;
run;



proc sort data=sug_all_ej_by_gene2; by gene_id; run;
proc sort data=sug_all_eu_by_gene2; by gene_id; run;
proc sort data=sug_all_ir_by_gene2; by gene_id; run;
proc sort data=sug_all_es_by_gene2; by gene_id; run;
proc sort data=sug_all_ad_by_gene2; by gene_id; run;
proc sort data=sug_all_aa_by_gene2; by gene_id; run;
proc sort data=sug_all_ada_by_gene2; by gene_id; run;
proc sort data=sug_dtct_ej_by_gene2; by gene_id; run;
proc sort data=sug_dtct_eu_by_gene2; by gene_id; run;
proc sort data=sug_dtct_ir_by_gene2; by gene_id; run;
proc sort data=sug_dtct_es_by_gene2; by gene_id; run;
proc sort data=sug_dtct_ad_by_gene2; by gene_id; run;
proc sort data=sug_dtct_aa_by_gene2; by gene_id; run;
proc sort data=sug_dtct_ada_by_gene2; by gene_id; run;
proc sort data=sug_tested_ej_by_gene2; by gene_id; run;
proc sort data=sug_tested_eu_by_gene2; by gene_id; run;
proc sort data=sug_tested_ir_by_gene2; by gene_id; run;
proc sort data=sug_tested_es_by_gene2; by gene_id; run;
proc sort data=sug_tested_ad_by_gene2; by gene_id; run;
proc sort data=sug_tested_aa_by_gene2; by gene_id; run;
proc sort data=sug_tested_ada_by_gene2; by gene_id; run;
proc sort data=sug_de_ej_by_gene2; by gene_id; run;
proc sort data=sug_de_eu_by_gene2; by gene_id; run;
proc sort data=sug_de_ir_by_gene2; by gene_id; run;
proc sort data=sug_de_es_by_gene2; by gene_id; run;
proc sort data=sug_de_ad_by_gene2; by gene_id; run;
proc sort data=sug_de_aa_by_gene2; by gene_id; run;
proc sort data=sug_de_ada_by_gene2; by gene_id; run;


data splicing_by_gene_summary;
merge sug_all_ej_by_gene2 sug_all_eu_by_gene2 sug_all_ir_by_gene2 sug_all_es_by_gene2 sug_all_ad_by_gene2 sug_all_aa_by_gene2 sug_all_ada_by_gene2 sug_dtct_ej_by_gene2 sug_dtct_eu_by_gene2 sug_dtct_ir_by_gene2 sug_dtct_es_by_gene2 sug_dtct_ad_by_gene2 sug_dtct_aa_by_gene2 sug_dtct_ada_by_gene2 sug_tested_ej_by_gene2 sug_tested_eu_by_gene2 sug_tested_ir_by_gene2 sug_tested_es_by_gene2 sug_tested_ad_by_gene2 sug_tested_aa_by_gene2 sug_tested_ada_by_gene2 sug_de_ej_by_gene2 sug_de_eu_by_gene2 sug_de_ir_by_gene2 sug_de_es_by_gene2 sug_de_ad_by_gene2 sug_de_aa_by_gene2 sug_de_ada_by_gene2 ;
by gene_id;
run;


data sugrue.splicing_by_gene_summary2;
   set splicing_by_gene_summary;
   array change _numeric_;
            do over change;
            if change=. then change=0;
            end;
   run ;




