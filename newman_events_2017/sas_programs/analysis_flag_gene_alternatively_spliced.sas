/* MISO analysis:
1. Check that the set of genes with MISO SE are the same set with SE in Event analysis
2. Check that the set of genes with detected MISO SE are the same set with detected SE in Event analysis
3. Check that the set of genes with DE MISO SE are the same set with DE/DD SE in Event analysis */

ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Count the genes that have fragments/exons and/or junctions that are differentially
   detected between OLDs and NPCs. This is evidence of AS */

data frags_on;
  set event.flag_fragment_on;
  if sum(flag_fragment_nsc_on,flag_fragment_old_on) = 1;
  keep fragment_id;
run;

data events_on;
  set event.flag_splicing_on;
  if sum(flag_event_nsc_on,flag_event_old_on) = 1;
  keep event_id;
run;

data fus_on;
  set event.flag_fusion_on;
  if sum(flag_fusion_nsc_on,flag_fusion_old_on) = 1;
  keep fusion_id;
run;

data frag2gene;
  set mm10.mm10_fragment2exon2gene;
  keep fragment_id gene_id;
run;

data fus2gene;
  set mm10.mm10_refseq_fusion_si_info_v2;
  keep fusion_id primary_gene_id;
  rename primary_gene_id=gene_id;
run;

data event2gene;
  set evspl.splicing_events_annot_refseq;
  keep event_id gene_id;
run;

proc sort data=frags_on;
  by fragment_id;
proc sort data=events_on;
  by event_id;
proc sort data=fus_on;
  by fusion_id;
proc sort data=frag2gene nodup;
  by fragment_id gene_id;
proc sort data=fus2gene nodup;
  by fusion_id gene_id;
proc sort data=event2gene nodup;
  by event_id gene_id;
run;

data genes_w_dd_frag;
  merge frag2gene (in=in1) frags_on (in=in2);
  by fragment_id;
  if in1 and in2;
run;

data genes_w_dd_fus;
  merge fus2gene (in=in1) fus_on (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data genes_w_dd_event;
  merge event2gene (in=in1) events_on (in=in2);
  by event_id;
  if in1 and in2;
run;


data genes_w_dd_frag2;
  set genes_w_dd_frag;
  keep gene_id;
run;

data genes_w_dd_fus2;
  set genes_w_dd_fus;
  keep gene_id;
run;

data genes_w_dd_event2;
  set genes_w_dd_event;
  keep gene_id;
run;

data genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

proc sort data=genes nodup;
   by gene_id;
proc sort data=genes_w_dd_frag2 nodup;
   by gene_id;
proc sort data=genes_w_dd_fus2 nodup;
   by gene_id;
proc sort data=genes_w_dd_event2 nodup;
   by gene_id;
run;

data flag_gene_as;
  merge genes (in=in1) genes_W_dd_frag2 (in=in2) genes_W_dd_fus2 (in=in3) genes_W_dd_event2 (in=in4) ;
  by gene_id;
  if in2 then flag_gene_fragment_dd=1; else flag_Gene_Fragment_dd=0;
  if in3 then flag_gene_fusion_dd=1; else flag_Gene_fusion_dd=0;
  if in4 then flag_gene_splicing_dd=1; else flag_Gene_splicing_dd=0;
  if in1 then output;
run;

proc freq data=flag_gene_as noprint;
   tables flaG_gene_fragment_dd*flag_Gene_fusion_dd*flag_gene_splicing_dd / out=as_count;
run;

proc print data=as_count;
run;

/*
flag_gene_                  flag_gene_
 fragment_    flag_gene_     splicing_
    dd         fusion_dd        dd        COUNT

     0             0             0         3529
     0             0             1         4774
     1             0             0          303
     1             0             1         1416
     1             1             0         4678
     1             1             1         6175

*/

data flag_gene_as2;
  set flag_gene_as;
  if sum(flag_gene_fragment_dd,flag_gene_fusion_dd,flag_gene_splicing_dd) > 0
     then flag_gene_as=1;
     else flag_gene_as=0;
run;

/* Make permenant */

data event.flag_gene_alt_spliced;
  set flag_gene_as2;
run;

proc freq data=event.flag_gene_alt_spliced;
   tables flag_Gene_as;
run;


/*
                                           Cumulative    Cumulative
  flag_gene_as    Frequency     Percent     Frequency      Percent
  -----------------------------------------------------------------
             0        3529       16.91          3529        16.91
             1       17346       83.09         20875       100.00


*/

