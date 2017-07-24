ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that the set of genes with detected MISO SE are the same set with detected SE in Event analysis
  For this, I only want to look at genes that both a MISO event and an exonskipping junction */

data miso_se;
   set event.miso_events_gene_in_refseq;
   keep ens_gene_id event_name;
run;

data miso_analyzed;
   set event.miso_all_results_nsc_v_old_trim;
   where NSC_OLD_diff ne .;
   keep event_name;
run;

proc sort data=miso_se;
  by event_name;
proc sort data=miso_analyzed;
  by event_name;
run;

data miso_analyzed2ens;
  merge miso_se (in=in1) miso_analyzed (in=in2);
  by event_name;
  if in1 and in2;
run; *2530 events analyzed;

data miso_genes;
  set miso_analyzed2ens;
  keep ens_gene_id;
run;

proc sort data=miso_genes nodup;
  by ens_gene_id;
run;

data ens2refseq;
   set event.ensembl2refseq_gene_id;
   keep ens_gene_id gene_id;
run;

proc sort data=miso_genes nodup;
   by ens_gene_id;
proc sort data=ens2refseq nodup;
   by ens_gene_id gene_id;
run;

data miso_se2refseq_dtct;
   merge miso_genes (in=in1) ens2refseq (in=in2);
   by ens_gene_id;
   if in1 and in2 then flag_has_refseq=1;
   else if in1 then flag_has_refseq=0;
   if in1 then output;
run;  *1667 RefSeq genes;


*Exonskipping genes;

data exonskip;
  set evspl.splicing_events_annot_refseq;
  where flag_exonskip=1;
  keep gene_id event_id;
run;

data event_dtct;
   set event.flag_splicing_on;
   where flag_event_nsc_on=1 or flag_event_old_on=1;
run;

proc sort data=exonskip;
   by event_id;
proc sort data=event_dtct;
   by event_id;
run;

data exonskip_dtct;
   merge exonskip (in=in1) event_dtct (in=in2);
   by event_id;
   if flag_event_nsc_on=1 and flag_event_old_on=1 then flag_event_on_both=1;
   else flag_event_on_both=0;
   flag_event=1;
run;

data exonskip_dtct_both exonskip_dtct_one;
   set exonskip_dtct;
   if flag_Event_on_both=1 then output exonskip_dtct_both;
   else output exonskip_dtct_one;
   keep gene_id;
run;

proc sort data=exonskip_dtct_both nodup;
   by gene_id;
proc sort data=exonskip_dtct_one nodup;
   by gene_id;
proc sort data=miso_se2refseq_dtct;
   by gene_id;
run;

data miso2es_dtct;
  merge miso_se2refseq_dtct (in=in1) exonskip_dtct_both (in=in2) exonskip_dtct_one (in=in3);
   by gene_id;
  if in1 then flag_has_miso_se_dtct=1; else flag_has_miso_se_dtct=0;
  if in2 then flag_exonskip_dtct_both=1; else flag_exonskip_dtct_both=0;
  if in3 then flag_exonskip_dtct_one=1; else flag_exonskip_dtct_one=0;
run;

*MISO genes;
data genes_to_keep;
   set event.miso_refseq_exonskip_Compare;
   where flag_has_refseq=1 and flag_has_miso_se=1 and flag_has_exonskip_event=1;
   keep gene_id;
run; *2014 genes;

proc sort data=genes_to_keep nodup;
   by gene_id;
proc sort data=miso2es_dtct nodup;
   by gene_id;
run;

data miso2es_dtct2;
  merge miso2es_dtct (in=in1) genes_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
run;




proc freq data=miso2es_dtct2 noprint;
  tables flag_has_miso_se_dtct*flag_exonskip_dtct_both*flag_exonskip_dtct_one / out=es_count;
run;

proc print data=es_count;
run;

/*
flag_has_
 miso_se_    flag_exonskip_    flag_exonskip_
   dtct         dtct_both         dtct_one       COUNT

    0               0                 1           283
    0               1                 1           139
    1               0                 1           738
    1               1                 0             7
    1               1                 1           846


283 genes with ES junctions only seen in either NPCs and OLDs
139 genes with ES junctions in NPCS and OLDs and ES junctions only seen in either NPCs and OLDs
738 genes with MISO events analyzed and ES junctions only seen in either NPCs and OLDs
7 genes with MISO events analyzed and ES junctions in NPCS and OLDs
846 in all

*/
/* Make permenant */

data event.miso_refseq_exonskip_cmpr_dtct;
   set miso2es_dtct2;
run;

 

